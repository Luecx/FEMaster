/**
* @file command_spec.cpp
 * @brief Implements the CommandSpec struct used for parser command registration.
 *
 * @date 06.10.2025
 * @version 1.0
 */

#include "command_spec.h"

namespace fem::reader2 {

CommandSpec& CommandSpec::doc(std::string summary, std::string notes) {
    doc_summary = std::move(summary);
    doc_notes   = std::move(notes);
    return *this;
}

CommandSpec& CommandSpec::plan(Condition cnd,
                               std::initializer_list<Segment> segs,
                               std::string display_name) {
    ConditionalPlan p = ConditionalPlan::make()
                            .when(std::move(cnd))
                            .segments(segs)
                            .name(std::move(display_name));
    plans.push_back(std::move(p));
    return *this;
}

CommandSpec& CommandSpec::plan(Condition cnd,
                               Segment segment,
                               std::string display_name) {
    ConditionalPlan p = ConditionalPlan::make()
                            .when(std::move(cnd))
                            .segments({segment})
                            .name(std::move(display_name));
    plans.push_back(std::move(p));
    return *this;
}

CommandSpec& CommandSpec::opens_scope(Scope s) {
    open_scope = std::move(s);
    return *this;
}

CommandSpec& CommandSpec::keys(KeyRules k) {
    key_rules = std::move(k);
    return *this;
}

CommandSpec& CommandSpec::on_enter(std::function<void(Context&, const Keyword&)> f) {
    on_enter_cb = std::move(f);
    return *this;
}

CommandSpec& CommandSpec::on_keyword(std::function<void(Context&, const Keyword&)> f) {
    on_keyword_cb = std::move(f);
    return *this;
}

CommandSpec& CommandSpec::on_exit(std::function<void(Context&)> f) {
    on_exit_cb = std::move(f);
    return *this;
}

} // namespace fem::reader2
