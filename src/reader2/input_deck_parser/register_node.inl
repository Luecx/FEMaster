// commands/register_node.inl
// Intentionally no include guards — include this in exactly ONE translation unit.

namespace fem::reader2::commands {

inline void register_node(Registry& reg, fem::model::Model& model, reader::Writer& writer) {
    (void)writer; // reserved for future use

    reg.register_command("ROOT", "NODE",
        [&model](CommandSpec& c) {
            // Keys: optional NSET with default "NALL"
            c.keys(
                KeyRules::make()
                    .optional("NSET", "NALL")
                    .desc("Active node set")
            );

            // Runs once per *NODE keyword line (before any data segments).
            c.on_enter([&model](Context&, const Keyword& kw) {
                model._data->node_sets.activate(kw.get<std::string>("NSET"));
            });

            // Data lines: ID, then optional X, Y, Z (each defaulting to 0 if omitted/empty).
            c.pattern(
                Pattern::make()
                    .fixed<int, 1>("ID").desc("node id")

                    .fixed<Precision, 1>("X")
                        .allow_missing(true).on_na_value(0).on_empty_value(0)
                        .desc("x-coordinate (optional, default 0)")

                    .fixed<Precision, 1>("Y")
                        .allow_missing(true).on_na_value(0).on_empty_value(0)
                        .desc("y-coordinate (optional, default 0)")

                    .fixed<Precision, 1>("Z")
                        .allow_missing(true).on_na_value(0).on_empty_value(0)
                        .desc("z-coordinate (optional, default 0)")

                    .bind_kv<int, Precision, Precision, Precision>(
                        [&model](Context&, const Keyword::Map&, int id, Precision x, Precision y, Precision z, const LineMeta&) {
                            model.set_node(id, x, y, z);
                        })
            );
        }
    );
}

} // namespace fem::reader2::commands
