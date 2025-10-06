/**
 * @file pattern.cpp
 * @brief Implements non-template members of Pattern.
 */

#include "pattern.h"

namespace fem::reader2 {

Pattern Pattern::make() {
    return Pattern{};
}

Pattern& Pattern::allow_multiline(bool on) {
    _multiline = on;
    return *this;
}

Pattern& Pattern::desc(std::string description) {
    if (_elements.empty()) {
        throw std::runtime_error("Pattern::desc(): no element to describe (call fixed()/var() first)");
    }
    _elements.back().set_description(std::move(description));
    return *this;
}

Pattern& Pattern::summary(std::string text) {
    _summary = std::move(text);
    return *this;
}

Pattern& Pattern::notes(std::string text) {
    _notes = std::move(text);
    return *this;
}

bool Pattern::multiline() const {
    return _multiline;
}

void Pattern::validate() const {
    size_t var_seen = 0;
    for (size_t index = 0; index < _elements.size(); ++index) {
        if (_elements[index].type() == PatternElement::Type::Variable) {
            ++var_seen;
            if (index + 1 != _elements.size()) {
                throw std::runtime_error("Pattern: variable element must be last");
            }
        }
    }
    if (var_seen > 1) {
        throw std::runtime_error("Pattern: only one variable element supported");
    }
}

size_t Pattern::min_required_tokens() const {
    size_t count = 0;
    for (const auto& element : _elements) {
        count += element.min_count();
    }
    return count;
}

size_t Pattern::max_allowed_tokens() const {
    size_t count = 0;
    for (const auto& element : _elements) {
        count += element.max_count();
    }
    return count;
}

std::vector<std::string> Pattern::column_names() const {
    std::vector<std::string> names;
    auto push_sequence = [&](const std::string& base, size_t total) {
        if (total == 1) {
            names.push_back(base);
            return;
        }
        for (size_t i = 1; i <= total; ++i) {
            names.push_back(base + std::to_string(i));
        }
    };
    for (const auto& element : _elements) {
        push_sequence(element.base_name(), element.max_count());
    }
    return names;
}

const std::string& Pattern::doc_summary() const {
    return _summary;
}

const std::string& Pattern::doc_notes() const {
    return _notes;
}

const std::vector<PatternElement>& Pattern::elements() const {
    return _elements;
}

std::shared_ptr<void> Pattern::convert(const std::vector<std::string>& tokens,
                                       const std::vector<size_t>& counts) const {
    if (!_converter) {
        throw std::runtime_error("Pattern not bound");
    }
    return _converter(tokens, counts);
}

void Pattern::invoke(Context& context, const Keyword& kw, void* tuple_ptr, const LineMeta& meta) const {
    if (!_invoker_with_kv) {
        throw std::runtime_error("Pattern not bound (bind_kv missing)");
    }
    _invoker_with_kv(context, kw.kv(), tuple_ptr, meta);
}

size_t Pattern::arity() const {
    return _arity;
}

} // namespace fem::reader2
