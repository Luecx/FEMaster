#include "reader.h"
#include "registry.h"
#include "types.h"
#include "keyword.h"

#include <stdexcept>
#include <sstream>

namespace fem::reader2 {

//------------------------------------------------------------------------------
// basic line buffering
//------------------------------------------------------------------------------
Line& Reader::peek_effective(File& file) {
    if (!_look) {
        _look = file.next_effective();
    }
    return *_look;
}

bool Reader::ensure_data_line(File& file){
    Line& P = peek_effective(file);
    if (P.type() != LineType::DATA) return false;
    return true;
}

size_t Reader::tokens_available() const {
    if (!_look || _look->type()!=LineType::DATA) return 0;
    const auto& values = _look->values();
    if (_tokenIndex >= values.size()) return 0;
    return values.size() - _tokenIndex;
}

size_t Reader::take_tokens(size_t count,
                           std::vector<std::string>& out,
                           [[maybe_unused]] const std::string& file,
                           [[maybe_unused]] int& last_line_no)
{
    if (!_look || _look->type() != LineType::DATA) return 0;
    const auto& values = _look->values();
    size_t avail = tokens_available();
    size_t take = (count < avail ? count : avail);
    for (size_t index = 0; index < take; ++index) {
        out.push_back(values[_tokenIndex++]);
    }
    return take;
}

void Reader::maybe_pop_line(){
    if (!_look || _look->type()!=LineType::DATA) return;
    const auto& values = _look->values();
    if (_tokenIndex >= values.size()) {
        _look.reset();
        _tokenIndex = 0;
    }
}

//------------------------------------------------------------------------------
// ensure_scope_accepts_or_bubble()
//------------------------------------------------------------------------------
size_t Reader::ensure_scope_accepts_or_bubble(const std::string& name) {
    while (true) {
        const Scope& scope = _context.current_scope();

        if (Registry::instance().find(scope, name)) {
            return _context.depth() - 1;
        }

        if (scope == "ROOT") {
            throw std::runtime_error(
                "Unexpected command '" + name +
                "' at ROOT (no matching handler registered).");
        }

        if (!_blocks.empty() && _blocks.back().scope_name == scope) {
            auto blk = _blocks.back();
            _blocks.pop_back();
            if (blk.on_exit_cb) {
                blk.on_exit_cb(_context);
            }
        }

        _context.pop_scope();
    }
}

//------------------------------------------------------------------------------
// run()
//------------------------------------------------------------------------------
void Reader::run(const std::string& filepath) {
    File file(filepath);

    while (true) {
        Line& L = peek_effective(file);
        if (L.type() == LineType::EOF_MARK) break;

        if (L.type() != LineType::KEYWORD) {
            // ignore stray data/comments until a keyword
            _look.reset();
            _tokenIndex = 0;
            continue;
        }

        std::string name = L.command();

        // Close scopes until a scope accepts this command
        ensure_scope_accepts_or_bubble(name);

        const Scope& scope = _context.current_scope();
        const CommandSpec* spec = Registry::instance().find(scope, name);
        if (!spec) {
            _look.reset();
            _tokenIndex = 0;
            continue;
        }

        // Build validated Keyword
        Keyword kw = Keyword::from_line(L, scope, spec->key_rules);

        // consume keyword
        _look.reset();
        _tokenIndex = 0;

        // callbacks
        if (spec->on_enter_cb) {
            spec->on_enter_cb(_context, kw);
        }
        if (spec->on_keyword_cb) {
            spec->on_keyword_cb(_context, kw);
        }

        // pick plan
        const ConditionalPlan* plan = nullptr;
        for (const auto& p : spec->plans) {
            if (p.cond()(kw)) { plan = &p; break; }
        }
        if (!plan) {
            throw std::runtime_error("No plan matched for command '" + name + "' at scope '" + scope + "'");
        }

        // execute segments in order
        for (const auto& seg : plan->segments()) {
            const Pattern& pat = seg.pattern();
            pat.validate();

            size_t groups = 0;
            while (groups < seg.range().max()) {
                // Can we start a new group?
                Line& P = peek_effective(file);
                if (P.type() == LineType::EOF_MARK || P.type() == LineType::KEYWORD) break;
                if (P.type() != LineType::DATA) {
                    _look.reset();
                    _tokenIndex = 0;
                    break;
                }

                // gather tokens
                const size_t min_need = pat.min_required_tokens();
                const size_t max_allow = pat.max_allowed_tokens();

                std::vector<std::string> toks;
                toks.reserve(max_allow);
                LineMeta last_meta{file.path(), file.line_number()};
                int last_line_no = file.line_number();

                // Step 1: fulfill minimum
                size_t have = 0;
                if (!pat.multiline()) {
                    const auto& vals = P.values();
                    (void)vals;
                    size_t avail = tokens_available();
                    if (avail < min_need) {
                        if (groups >= seg.range().min()) break;
                        std::ostringstream oss;
                        oss << "Too few columns for pattern (need " << min_need
                            << ", have " << avail << ") at " << file.path() << ":" << file.line_number();
                        throw std::runtime_error(oss.str());
                    }
                    size_t t = take_tokens(min_need, toks, file.path(), last_line_no);
                    last_meta = {file.path(), static_cast<int>(last_line_no)};
                    have += t;
                } else {
                    // multiline: slurp across lines until min_need reached
                    while (have < min_need) {
                        if (!ensure_data_line(file)) break;
                        size_t t = take_tokens(min_need - have, toks, file.path(), last_line_no);
                        if (t == 0) break;
                        have += t;
                        last_meta = {file.path(), static_cast<int>(last_line_no)};
                        if (have < min_need) {
                            // need a new line
                            maybe_pop_line();
                            Line& nxt = peek_effective(file);
                            if (nxt.type() != LineType::DATA) break;
                        }
                    }
                    if (have < min_need) {
                        if (groups >= seg.range().min()) break;
                        std::ostringstream oss;
                        oss << "Pattern requires at least " << min_need << " tokens, got " << have
                            << " at " << file.path() << ":" << file.line_number();
                        throw std::runtime_error(oss.str());
                    }
                }

                // Step 2: optional extras only from the same last line
                size_t extras_cap = max_allow - min_need;
                size_t extras_take = 0;
                if (extras_cap > 0) {
                    extras_take = std::min(extras_cap, tokens_available());
                    if (extras_take) {
                        size_t t = take_tokens(extras_take, toks, file.path(), last_line_no);
                        (void)t;
                        last_meta = {file.path(), static_cast<int>(last_line_no)};
                    }
                }

                maybe_pop_line();

                // counts per element
                std::vector<size_t> counts;
                counts.reserve(pat.elements().size());
                for (const auto& element : pat.elements()) {
                    if (element.type() == PatternElement::Type::Variable) {
                        counts.push_back(element.effective_count(extras_take));
                    } else {
                        counts.push_back(element.effective_count(0));
                    }
                }

                auto tup_ptr = pat.convert(toks, counts);
                pat.invoke(_context, kw, tup_ptr.get(), last_meta);

                ++groups;
            }

            if (groups < seg.range().min()) {
                std::ostringstream oss;
                oss << "Too few groups for segment '" << (seg.label().empty()?std::string("Segment"):seg.label())
                    << "' (min=" << seg.range().min() << ")";
                throw std::runtime_error(oss.str());
            }
        }

        // scope mgmt (unchanged)
        const Scope open_as = spec->open_scope ? *spec->open_scope : spec->name;
        bool explicit_scope = static_cast<bool>(spec->open_scope);
        bool has_children   = Registry::instance().scope_has_children(open_as);

        if (explicit_scope) {
            _context.push_scope(open_as);
            _blocks.push_back(ActiveBlock{open_as, spec->on_exit_cb});
        } else if (has_children) {
            _context.push_scope(open_as);
            _blocks.push_back(ActiveBlock{open_as, spec->on_exit_cb});
        }
    }

    // close all blocks at EOF
    while (!_blocks.empty()) {
        auto blk = _blocks.back();
        _blocks.pop_back();
        if (blk.on_exit_cb) {
            blk.on_exit_cb(_context);
        }
        _context.pop_scope();
    }
}

} // namespace fem::reader2
