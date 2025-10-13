#pragma once
/**
 * @file segments.h
 * @brief Ein Block aus einem oder mehreren Segmenten, optional bewacht durch eine KeyCondition,
 *        die gegen das *aktuelle Keyword* evaluiert wird.
 */

#include "key_condition.h"
#include "keyword.h"
#include "segment.h"
#include "types.h"

#include <optional>
#include <vector>

namespace fem::reader2 {

/**
 * @brief Container für mehrere Segmente, optional mit Schlüsselbedingung.
 *
 * Wenn `when()` nicht gesetzt ist, gilt der Block immer.
 * Ist `when()` gesetzt, muss die Bedingung gegen `kw.kv()` wahr sein, sonst wird der Block übersprungen.
 */
class Segments {
    public:
    /// Fabrik
    static Segments make() { return Segments{}; }

    /// Bedingung für diesen Block (gegen das aktuelle Keyword geprüft)
    Segments& when(condition::KeyCondition kc) {
        _when = std::move(kc);
        return *this;
    }

    /// Ein Segment hinzufügen
    Segments& add(Segment s) {
        _items.push_back(std::move(s));
        return *this;
    }

    /// Mehrere Segmente hinzufügen
    Segments& add(std::initializer_list<Segment> list) {
        _items.insert(_items.end(), list.begin(), list.end());
        return *this;
    }

    /// True, wenn der Block für das gegebene Keyword gilt (oder keine Bedingung gesetzt ist).
    bool matches(const Keyword& kw) const {
        if (!_when) return true;
        return _when->eval(kw.kv());
    }

    /// Zugriff
    const std::vector<Segment>& items() const { return _items; }
    bool has_condition() const { return _when.has_value(); }

    private:
    std::vector<Segment>                   _items;
    std::optional<condition::KeyCondition> _when;
};

} // namespace fem::reader2
