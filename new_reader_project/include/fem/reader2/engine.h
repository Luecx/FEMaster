#pragma once
#include <string>
#include <vector>
#include <stdexcept>
#include <sstream>
#include <iostream>

#include "registry.h"
#include "../reader/file.h"
#include "../reader/line.h"

namespace fem { namespace reader2 {

inline Keys keys_from_line(const fem::reader::Line& l){
    Keys k;
    // Normalisierte Keyword-Zeile: COMMAND,KEY=VAL,FLAG
    std::string s = l.line();
    std::stringstream ss(s);
    std::string cmd;
    std::getline(ss, cmd, ',');
    std::string item;
    while (std::getline(ss, item, ',')){
        auto pos = item.find("=");
        if (pos==std::string::npos) k.kv[item] = "";
        else k.kv[item.substr(0,pos)] = item.substr(pos+1);
    }
    return k;
}

inline std::size_t pattern_required_tokens(const Pattern& p){
    std::size_t n=0; for (auto& e : p.elems) n += e->count(); return n;
}

inline void append_values_upper(const fem::reader::Line& l, std::vector<std::string>& out){
    for (auto& v : l.values()) out.push_back(v);
}

struct Engine {
    const Registry& reg;
    explicit Engine(const Registry& r): reg(r) {}

    void run(fem::reader::File& file){
        using LT = fem::reader::LineType;

        // Scope-Stack: nur der direkte Parent ist relevant, aber wir können „hochsteigen“
        std::vector<ParentInfo> scope;
        scope.push_back( ParentInfo{ "ROOT", Keys{} } );

        while (true){
            auto& ln = file.next_line();
            if (ln.type()==LT::END_OF_FILE) break;
            if (ln.type()==LT::DATA_LINE) {
                throw std::runtime_error("Unexpected DATA line without active command");
            }
            if (ln.type()!=LT::KEYWORD_LINE) continue;

            std::string cmd = ln.command();
            Keys self_keys  = keys_from_line(ln);

            // Command registriert?
            auto* spec = reg.find(cmd);

            // Wenn nicht registriert: Wir behandeln es als strukturierendes Parent-Keyword.
            // -> push als neuer Parent (der direkte Parent ist die aktuelle Scope-Spitze).
            if (!spec){
                scope.push_back( ParentInfo{cmd, self_keys} );
                continue;
            }

            // Registriert: Prüfe Zulassung. Falls nicht zugelassen, steige im Scope hoch.
            // Wir suchen den ersten (nächsten) Vorfahren (einschließlich aktueller Spitze),
            // unter dem die admit_-Bedingung erfüllt ist (oder keine admit_ vorhanden ist).
            int chosen_parent_index = -1;
            for (int i = static_cast<int>(scope.size()) - 1; i >= 0; --i) {
                const ParentInfo& candidate = scope[static_cast<std::size_t>(i)];
                if (!spec->admit_ || spec->admit_->eval(candidate, self_keys)) {
                    chosen_parent_index = i;
                    break;
                }
            }

            if (chosen_parent_index < 0){
                // Kein zulässiger Vorfahr → harter Fehler mit Diagnose
                std::ostringstream os;
                os << "Command '" << cmd << "' not admitted in current scope (";
                os << "stack: ";
                for (std::size_t i=0;i<scope.size();++i){
                    if (i) os << " > ";
                    os << scope[i].command;
                }
                os << ")";
                throw std::runtime_error(os.str());
            }

            // Scope auf die gefundene Ebene stutzen (hochgestiegen)
            scope.resize(static_cast<std::size_t>(chosen_parent_index + 1));

            // Variante wählen
            const Variant* chosen = nullptr;
            for (auto& v : spec->variants_){
                if (!v.condition_ || v.condition_->eval(scope.back(), self_keys)) {
                    chosen = &v;
                    break;
                }
            }
            if (!chosen) throw std::runtime_error("No matching variant for command: " + cmd);

            // Segmente ausführen
            for (auto& seg : chosen->segments_){
                std::vector<std::string> tokens;

                std::size_t lines_needed_min  = seg.range_.min_;
                std::size_t lines_allowed_max = seg.range_.max_;
                std::size_t lines_read        = 0;
                std::size_t needed_tokens     = pattern_required_tokens(seg.pattern_);

                while (true){
                    auto& dl = file.next();
                    if (dl.type()==LT::END_OF_FILE)
                        throw std::runtime_error("Unexpected EOF while reading segment for " + cmd);
                    if (dl.type()==LT::KEYWORD_LINE)
                        throw std::runtime_error("Encountered next command while reading data segment for " + cmd);
                    if (dl.type()==LT::DATA_LINE){
                        append_values_upper(dl, tokens);
                        ++lines_read;
                        if (!seg.pattern_.multiline) {
                            if (lines_read >= lines_needed_min) break;
                            if (lines_read >= lines_allowed_max) break;
                        } else {
                            if (tokens.size() >= needed_tokens) break;
                            if (lines_read >= lines_allowed_max) break;
                        }
                    }
                }

                if (tokens.size() < needed_tokens){
                    std::ostringstream os; os << "Insufficient tokens for segment of " << cmd
                       << " ("<<tokens.size()<<"/"<<needed_tokens<<")";
                    throw std::runtime_error(os.str());
                }
                if (tokens.size() > needed_tokens) tokens.resize(needed_tokens);

                if (!seg.invoke_) throw std::runtime_error("No binder for a segment of "+cmd);
                seg.invoke_(tokens);
            }

            // Dieses (erfolgreich verarbeitete) Keyword wird jetzt neuer Parent (Scope nach unten)
            scope.push_back( ParentInfo{cmd, self_keys} );
        }
    }
};

}} // namespace fem::reader2
