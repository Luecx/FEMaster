#pragma once

#include "section_shell.h"

namespace fem {

struct ABDShellSection : ShellSection {
   using Ptr = std::shared_ptr<ABDShellSection>;

   StaticMatrix<6, 6> abd = StaticMatrix<6, 6>::Zero();
   StaticMatrix<2, 2> shear = StaticMatrix<2, 2>::Zero();

   StaticMatrix<6, 6> get_abd() override;
   StaticMatrix<2, 2> get_shear() override;

   void info() override;
   std::string str() const override;
};

} // namespace fem
