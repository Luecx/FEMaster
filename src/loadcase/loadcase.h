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

    bool report_constraints = false;

    LoadCase(ID id, reader::Writer* writer, model::Model* model) : m_id(id), m_writer(writer), m_model(model) {}

    virtual void run() = 0;
    inline ID id() const { return m_id; }

protected:
    void report_constraint_groups(const model::Model::ConstraintGroups& groups) const;

};


}
}
