#pragma once

#include "../model/model.h"
#include "../reader/writer.h"

namespace fem{
namespace loadcase{

// every load case has loads and supports
struct LoadCase{

    const ID m_id;
    reader::Writer* m_writer;
    model::Model  * m_model;

    LoadCase(ID id, reader::Writer* writer, model::Model* model) : m_id(id), m_writer(writer), m_model(model) {}

    virtual void run() = 0;

};


}
}