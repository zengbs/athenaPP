//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file sr_blast.cpp
//  \brief Problem generator for spherical blast wave in flat spacetime.

// C headers

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // sqrt()
#include <cstring>    // strcmp()

// Athena++ headers
#include "../athena.hpp"                   // macros, enums
#include "../athena_arrays.hpp"            // AthenaArray
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"          // ParameterInput


Real threshold = 0.01;
 
int RefinementCondition(MeshBlock *pmb);
 
void Mesh::InitUserMeshData(ParameterInput *pin) {
  if (adaptive) {
    EnrollUserRefinementCondition(RefinementCondition);
  }
  return;
}



void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Prepare index bounds
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  if (block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  int kl = ks;
  int ku = ke;
  if (block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  // Get ratio of specific heats
  // Real gamma_adi = peos->GetGamma();
  // Real gamma_adi_red = gamma_adi / (gamma_adi - 1.0);

  // Read problem parameters
  Real center_x   = pin->GetReal("problem", "center_x"  );
  Real center_y   = pin->GetReal("problem", "center_y"  );
  Real center_z   = pin->GetReal("problem", "center_z"  );
  Real radius     = pin->GetReal("problem", "radius"    );
  Real rho_inner  = pin->GetReal("problem", "rho_inner" );
  Real tgas_inner = pin->GetReal("problem", "tgas_inner");
  Real rho_outer  = pin->GetReal("problem", "rho_outer" );
  Real tgas_outer = pin->GetReal("problem", "tgas_outer");

  AthenaArray<Real> b;
  b.NewAthenaArray(3, ncells3, ncells2, ncells1);

  // Initialize hydro variables
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {
        Real x1 = pcoord->x1v(i);
        Real x2 = pcoord->x2v(j);
        Real x3 = pcoord->x3v(k);
        Real separation = sqrt( SQR(x1-center_x) + SQR(x2-center_y) + SQR(x3-center_z) );


        // Set pressure and density
        Real rho, tgas;
        if (separation < radius) {
          rho = rho_inner;
          tgas = tgas_inner;
        } else {
          rho = rho_outer;
          tgas = tgas_outer;
        }

        // Set velocity
        Real ut = 1.0;
        Real ux = 0.0;
        Real uy = 0.0;
        Real uz = 0.0;
        Real u0, u1, u2, u3;

        u0 = ut; 
        u1 = ux; 
        u2 = uy; 
        u3 = uz; 


        phydro->w(IDN,k,j,i) = phydro->w1(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = phydro->w1(IPR,k,j,i) = tgas*rho;
        phydro->w(IVX,k,j,i) = phydro->w1(IVX,k,j,i) = u1 / u0;
        phydro->w(IVY,k,j,i) = phydro->w1(IVY,k,j,i) = u2 / u0;
        phydro->w(IVZ,k,j,i) = phydro->w1(IVZ,k,j,i) = u3 / u0;

        b(IB1,k,j,i) = 0.0;
        b(IB2,k,j,i) = 0.0;
        b(IB3,k,j,i) = 0.0;
      }
    }
  }
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, il, iu, jl, ju, kl, ku);

  return;
}


int RefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps = 0.0;

// 3D flag
  if (pmb->pmy_mesh->f3) {
    for (int k=pmb->ks-1; k<=pmb->ke+1; k++) {
      for (int j=pmb->js-1; j<=pmb->je+1; j++) {
        for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
          Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                               +SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i)))
                               +SQR(0.5*(w(IPR,k+1,j,i) - w(IPR,k-1,j,i))))/w(IPR,k,j,i);
          maxeps = std::max(maxeps, eps);
        }
      }   
    }   
// 2D flag
  } else if (pmb->pmy_mesh->f2) {
    int k = pmb->ks;
    for (int j=pmb->js-1; j<=pmb->je+1; j++) {
      for (int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real eps = std::sqrt(SQR(0.5*(w(IPR,k,j,i+1) - w(IPR,k,j,i-1)))
                             + SQR(0.5*(w(IPR,k,j+1,i) - w(IPR,k,j-1,i))))/w(IPR,k,j,i);
        maxeps = std::max(maxeps, eps);
      }   
    }   
  } else {
    return 0;
  }
 
  if (maxeps > threshold) return 1;
  if (maxeps < 0.25*threshold) return -1; 
  return 0;                                                                                                                              
}

