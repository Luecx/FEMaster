// parser/keys.h
#pragma once
#include <string>
#include <unordered_map>
#include <sstream>
#include <stdexcept>
#include <type_traits>

namespace fem::reader2 {

    struct KeyBag {
        std::unordered_map<std::string,std::string> kv;

        bool has(const std::string& k) const { return kv.find(k)!=kv.end(); }

        template<class T>
        static T convert(const std::string& s){
            if constexpr (std::is_same_v<T, std::string>) {
                return s;
            } else {
                std::istringstream iss(s);
                T v{};
                iss >> v;
                if (iss.fail()) throw std::runtime_error("Key conversion failed for value '"+s+"'");
                return v;
            }
        }

        template<class T>
        T get(const std::string& k, const T& def) const {
            auto it = kv.find(k);
            if (it==kv.end()) return def;
            return convert<T>(it->second);
        }

        template<class T>
        T require(const std::string& k) const {
            auto it = kv.find(k);
            if (it==kv.end()) throw std::runtime_error("Missing required key: "+k);
            return convert<T>(it->second);
        }
    };

} // namespace fem::reader2
