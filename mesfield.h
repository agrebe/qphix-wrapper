#pragma once

#include "qdp.h"

namespace QPhiX
{
//! Calculates the antihermitian field strength tensor  iF(mu,nu)
/*
 * \ingroup glue
 *
 *    F(mu,nu) =  (1/4) sum_p (1/2) [ U_p(x) - U^dag_p(x) ]
 *
 *  where
 *    U_1 = u(x,mu)*u(x+mu,nu)*u_dag(x+nu,mu)*u_dag(x,nu)
 *    U_2 = u(x,nu)*u_dag(x-mu+nu,mu)*u_dag(x-mu,nu)*u(x-mu,mu)
 *    U_3 = u_dag(x-mu,mu)*u_dag(x-mu-nu,nu)*u(x-mu-nu,mu)*u(x-nu,nu)
 *    U_4 = u_dag(x-nu,nu)*u(x-nu,mu)*u(x-nu+mu,nu)*u_dag(x,mu)

 *  NOTE: I am using  i*F  of the usual F defined by UKQCD, Heatlie et.al.

 * Arguments:

 *  \param f   field strength tensor f(cb,mu,nu) (Write)
 *  \param u   gauge field (Read)
 */

void mesField(QDP::multi1d<QDP::LatticeColorMatrixF> &f,
              const QDP::multi1d<QDP::LatticeColorMatrixF> &u);

void mesField(QDP::multi1d<QDP::LatticeColorMatrixD> &f,
              const QDP::multi1d<QDP::LatticeColorMatrixD> &u);
}
