#pragma once
/**
 * @file key_condition.h
 * @brief Boolescher Ausdruck über Schlüsseltests (present / missing / eq / in / neq / nin)
 *        zur Auswertung gegen eine Keyword-Map (z. B. des aktuellen Keywords oder Scopes).
 */

#include "keyword.h"

#include <initializer_list>
#include <string>
#include <unordered_set>
#include <vector>

namespace fem::reader2::condition {

/**
 * @brief Elementarer Schlüsseltest (Blatt eines booleschen Ausdrucksbaums).
 */
struct KeyTest {
    enum class Type { Present, Missing, Eq, In, Neq, Nin };

    std::string                     key;    ///< Schlüsselname (uppercased wie im Parser)
    Type                            type;   ///< Testtyp
    std::unordered_set<std::string> values; ///< Vergleichswerte (für Eq/In/Neq/Nin)

    // Fabriken
    static KeyTest present(std::string k);
    static KeyTest missing(std::string k);
    static KeyTest eq(std::string k, std::string v);
    static KeyTest in(std::string k, std::initializer_list<std::string> vs);
    static KeyTest neq(std::string k, std::string v);
    static KeyTest nin(std::string k, std::initializer_list<std::string> vs);
};

/**
 * @brief Boolescher Ausdruck über @ref KeyTest (AND / OR / NOT).
 *        Evaluiert gegen eine einzelne Key-Map.
 */
class KeyCondition {
public:
    enum class Op { Leaf, And, Or, Not };

    // Fabriken (Struktur)
    static KeyCondition Leaf(KeyTest t);
    static KeyCondition And(std::initializer_list<KeyCondition> xs);
    static KeyCondition Or (std::initializer_list<KeyCondition> xs);
    static KeyCondition Not(KeyCondition x);

    // Kurzformen (Leaves)
    static KeyCondition present(std::string k) { return Leaf(KeyTest::present(std::move(k))); }
    static KeyCondition missing(std::string k) { return Leaf(KeyTest::missing(std::move(k))); }
    static KeyCondition eq(std::string k, std::string v) { return Leaf(KeyTest::eq(std::move(k), std::move(v))); }
    static KeyCondition in(std::string k, std::initializer_list<std::string> vs) { return Leaf(KeyTest::in(std::move(k), vs)); }
    static KeyCondition neq(std::string k, std::string v) { return Leaf(KeyTest::neq(std::move(k), std::move(v))); }
    static KeyCondition nin(std::string k, std::initializer_list<std::string> vs) { return Leaf(KeyTest::nin(std::move(k), vs)); }

    /// Evaluiert gegen eine Key-Map.
    bool eval(const fem::reader2::Keyword::Map& kv) const;

private:
    Op                              _op { Op::Leaf };
    KeyTest                         _test;
    std::vector<KeyCondition>       _children;
};

} // namespace fem::reader2::condition
