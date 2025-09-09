//
// Created by f_eggers on 04.09.2025.
//

#ifndef IP_DATA_DICT_H
#define IP_DATA_DICT_H



#include "data_storage.h"
#include "dict.h"

#include "../core/types_eig.h"
namespace fem::model {

enum IPDataEntries : Index {
    STRESS,
};

using IPDataDict = DataStorage<IPData>;
}
#endif //IP_DATA_DICT_H
