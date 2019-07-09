/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2018 Robin Trunk
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

/* youngLaplace3d.cpp
 * In this example a Young-Laplace test is performed. A spherical domain
 * of fluid 2 is immersed in fluid 1. A diffusive interface forms and the
 * surface tension can be calculated using the Laplace pressure relation.
 * The pressure difference is calculated between a point in the middle of
 * the circular domain and a point furthest away from it in the
 * computational domain (here left bottom corner).
 *
 * This example shows the simplest case for the free-energy model with two
 * fluid components.
 */

#include "olb3D.h"
#include "olb3D.hh"   // use only generic version!
#include <cstdlib>
#include <iostream>
#include <fstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;
#define DESCRIPTOR D3Q19<CHEM_POTENTIAL,FORCE>

// Parameters for the simulation setup
const T radius = 30;
const double N  = 0.25*radius;
const T nx = 12*radius;
const T ny = 8*radius;
const T nz = 4*radius;

const T alpha = 0.0379*radius; //1.5;     // Interfacial width         [lattice units]
const T kappages = (1/(0.0379*radius)) *(6/201.6);
const T kappa2 = kappages/2.5;//default 0.005;  // For surface tensions      [lattice units]
const T kappa1 = kappa2*1.5;//default 0.0075; // For surface tensions      [lattice units]
const T gama = 1.;       // For mobility of interface [lattice units]

//h_i neu (analog microFluidics2d) zur Einführung von freeEnergyBoundary
//h_i= Parameter related to resulting contact angle of the boundary. [lattice units]
const T h1 = 0.;                  // Contact angle 90 degrees   [lattice units]
const T h2 = 0.;                  // Contact angle 90 degrees   [lattice units]

const int maxIter  = 1000; //default 60.000 -> hier Simulationsschritte einstellen
const int vtkIter  = 200;
const int statIter = 200;


void prepareGeometry( SuperGeometry3D<T>& superGeometry,
                      UnitConverter<T,DESCRIPTOR> const& converter,
                      IndicatorF3D<T>& indicator)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // superGeometry.rename( 0, 2 ); --> DELETED, doesn't use indicator as passed argument from main.
  superGeometry.rename( 0, 2, indicator );
  superGeometry.rename( 2, 1, 0, 1, 1 );


  std::vector<T> origin( 3, T(0) );
  std::vector<T> extend( 3, T(0) );

  // BOUNDARY INDICATORS, alternativ direkt mit eps angeben
  T eps = converter.getPhysLength(1);   // PHYSICAL Delta_x
  T edge = 0.5*eps;                     //    WITHIN DOMAIN: INDICATOR BOX EDGE LIES IN BETWEEN NODES
  T safety = 2*eps;                     //    OUTSIDE OF DOMAIN: 2 NODES SAFETY OVERLAP

  // TOP WALL, origin y=ny oben, extend Ebene xz aufspannen
  origin[0] = 0.0 - safety;
  origin[1] = ny - edge;
  origin[2] = 0.0 - safety;
  extend[0] = nx + 2*safety;            // 2*safety, von -safe auf +safe
  extend[1] = edge + safety;
  extend[2] = nz + 2*safety;            // analog oben
  IndicatorCuboid3D<T> top( extend, origin );
  superGeometry.rename( 2, 3, top ); //ersetze MN=2 mit 3 für den indicator 'top'

  // BOTTOM WALL, orgin bleibt an Ort, extend Ebene xz aufspannen
  origin[1] = 0.0 - safety;
  IndicatorCuboid3D<T> bottom( extend, origin );
  superGeometry.rename( 2, 4, bottom ); //ersetze MN=2 mit 4 für den indicator 'top'

  // SIDE WALLS links+rechts
  origin[0] = 0.0 - safety;
  origin[1] = edge;
  origin[2] = 0.0 - safety;
  extend[0] = nx + 2*safety;
  extend[1] = ny - edge - eps;   // -eps, SINCE 2*edge LESS (DUE TO TOP/BOTTOM BOUNDARY IN CORNERS)
  extend[2] = edge + safety;
  IndicatorCuboid3D<T> wall1( extend, origin ); // left wall in x-direction
  superGeometry.rename( 2, 5, wall1 );
  //rechte Wand, origin um z=nz verschieben
  origin[2] = nz - edge;
  IndicatorCuboid3D<T> wall2( extend, origin ); // right wall in x-direction
  superGeometry.rename( 2, 6, wall2 );

  // NOTE: PERIODIC FRONT AND BACK WALL: MN = 1, BELONGS TO BULK.
  // Keine eigene MN für In/Outflow mehr, da in main periodic boundary

  // MATERIAL NUMBERS:
  //    BEFORE: INFLOW 3, OUTFLOW 4, TOP 5, BOTTOM 6
  //    NOW: BULK = 1, TOP = 3, BOTTOM = 4, SIDE WALL LEFT = 5, SIDE WALL RIGHT = 6
  //    [RENAMING CHECK = 2 (if not appearing in terminal output, every boundary node got hit by an indicator!)]


  // // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}






void prepareLattice( SuperLattice3D<T, DESCRIPTOR>& sLattice1,
                     SuperLattice3D<T, DESCRIPTOR>& sLattice2,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics1,
                     Dynamics<T, DESCRIPTOR>& bulkDynamics2,
                     UnitConverter<T, DESCRIPTOR>& converter,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBC1,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBC2,
                     SuperGeometry3D<T>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  //neu für omega
  T omega = converter.getLatticeRelaxationFrequency();

  //define lattice Dynamics, jeweils beide Lattices beachten
  //bulkDynamics1 = forced ForcedBGKdynamics
  //bulkDynamics2 = FreeEnergyBGKdynamics

  //MN=0 -> no dynamics
  sLattice1.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  //MN=1 Fluid hat bulkdynamics
  sLattice1.defineDynamics( superGeometry, 1, &bulkDynamics1 );
  sLattice2.defineDynamics( superGeometry, 1, &bulkDynamics2 );

  //MN=3 TOP
  sLattice1.defineDynamics( superGeometry, 3, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 3, &instances::getNoDynamics<T, DESCRIPTOR>() );

  //MN= 4 BOTTOM
  sLattice1.defineDynamics( superGeometry, 4, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 4, &instances::getNoDynamics<T, DESCRIPTOR>() );

  //MN=5 LEFT WALL
  sLattice1.defineDynamics( superGeometry, 5, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 5, &instances::getNoDynamics<T, DESCRIPTOR>() );

  //MN= 4 RIGHT WALL
  sLattice1.defineDynamics( superGeometry, 6, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 6, &instances::getNoDynamics<T, DESCRIPTOR>() );

  // STANDARD WALL FOR FREE ENERGY.
  //    NOTE: THIS IS A NO SLIP BOUNDARY (default) -> boundaryInstatiotor.h ändern
  //    LEAVE THAT FOR NOW SINCE MODEL IS CORRECT.
  //    MAYBE CHANGE LATER IF POSSIBLE, OTHERWISE MENTION THAT IN THESIS (CURRENTLY NO FREE ENERGY SLIP BOUNDARY IN OPENLB)

  //add free Energy boundaries für Wände
  //gemäß Doxygen für 3D, 2Phasen: x(superGeometry, MN, alpha, kappa1, kappa2, h1, h2, latticeNumber)
  sOnBC1.addFreeEnergyWallBoundary( superGeometry, 5, alpha, kappa1, kappa2, h1, h2, 1 );
  sOnBC2.addFreeEnergyWallBoundary( superGeometry, 5, alpha, kappa1, kappa2, h1, h2, 2 );

  sOnBC1.addFreeEnergyWallBoundary( superGeometry, 6, alpha, kappa1, kappa2, h1, h2, 1 );
  sOnBC2.addFreeEnergyWallBoundary( superGeometry, 6, alpha, kappa1, kappa2, h1, h2, 2 );


  // MOVING WALLS (top+bottom): use inlet velocity to fix velocity here in normal tangential direction
  auto topIndicator = superGeometry.getMaterialIndicator(3);
  sOnBC1.addFreeEnergyInletBoundary( topIndicator, omega, "velocity", 1 );
  sOnBC2.addFreeEnergyInletBoundary( topIndicator, omega, "velocity", 2 );

  auto bottomIndicator = superGeometry.getMaterialIndicator(4);
  sOnBC1.addFreeEnergyInletBoundary( bottomIndicator, omega, "velocity", 1 );
  sOnBC2.addFreeEnergyInletBoundary( bottomIndicator, omega, "velocity", 2 );

  //Geschwindigkeitsvektor
  std::vector<T> v( 3,T() );
  AnalyticalConst3D<T,T> zeroVelocity( v );   // NullGeschwindigkeit
  AnalyticalConst3D<T,T> zero ( 0. );         // Null
  AnalyticalConst3D<T,T> one ( 1. );          // Eins

  // TEST DIFFERENT COEFFICIENTS BEFORE ALPHA? (THIS MODIFIES THE INTERFACE THICKNESS!!)
  // T interfaceCoeff = 10.; // DEFAULT
  // Ziel:1
  //T interfaceCo = 5.;
  T interfaceCo = 10.;
  SmoothIndicatorSphere3D<T,T> sphere( {nx/2., ny/2., nz/2.}, radius, interfaceCo*alpha ); //Tropfen
  AnalyticalIdentity3D<T,T> rho( one ); //rho=1
  AnalyticalIdentity3D<T,T> phi( one - sphere - sphere ); //phi






  /*Davor:
  bei sLattice 1 mit rho, sLattice2 mit phi
  (bulkIndicator, rhoF, uF) ersetzt (superGeometry, 1, rho, zeroVelocity)
  auto bulkIndicator = superGeometry.getMaterialIndicator({ 1 }); //BULK IS SOLELY ON 1!
  sLattice1.iniEquilibrium( bulkIndicator, rho, zeroVelocity );
  sLattice2.iniEquilibrium( bulkIndicator, phi, zeroVelocity );
  */

  // OR: INITIALIZE ANYTHING WITH EQUILIBRIUM. (AND DEFINE RHO AND U ONLY ON MOVING WALLS)
  auto allIndicator = superGeometry.getMaterialIndicator({ 1, 3, 4, 5, 6 });
  sLattice1.iniEquilibrium( allIndicator, rho, zeroVelocity );
  sLattice2.iniEquilibrium( allIndicator, phi, zeroVelocity );

  //Geschwindigkeit für top, bottom wall
  clout << "getCharLatticeVelocity:" << converter.getCharLatticeVelocity() << std::endl;
  AnalyticalConst3D<T,T> uTop( converter.getCharLatticeVelocity(), T( 0 ), T( 0 ) );
  AnalyticalConst3D<T,T> uBottom( -converter.getCharLatticeVelocity(), T( 0 ), T( 0 ) );
  /*
  Alternative Überlegung ohne converter. Mit: uTop(v_x, v_y, v_z)
  v_x=1000 m/s entspräche:
  clout << T(1000.) << std::endl;
  AnalyticalConst3D<T,T> uTop( T(1000.) , T( 0 ), T( 0 ) );

  NOTE:
    charPhysVelocity is used in Re as u_w.
    converter.getCharLatticeVelocity() pulls charPhysVelocity and transforms it to lattice units
    if u_w gets changed, but charPhysVelocity not, then Re computation is false!
  */


  //MN=3 TOP
  sLattice1.defineU( topIndicator, uTop );
  sLattice2.defineU( topIndicator, uTop );
  sLattice1.defineRho( topIndicator, rho );
  sLattice2.defineRho( topIndicator, phi );

  //MN=4 BOTTOM
  sLattice1.defineU( bottomIndicator, uBottom );
  sLattice2.defineU( bottomIndicator, uBottom );
  sLattice1.defineRho( bottomIndicator, rho );
  sLattice2.defineRho( bottomIndicator, phi );

  // Make the lattice ready for simulation, initialise
  sLattice1.initialize();
  sLattice2.initialize();

  sLattice1.communicate();
  sLattice2.communicate();
  //--------------------------------------------------------------------------

  clout << "Prepare Lattice ... OK" << std::endl;
}

//------------------------------------------------------------------------------
//Hier war mal void setBoundaryValues -> jetzt weg bzw Inhalt oben in prepareLattice
//------------------------------------------------------------------------------


//Kopplung
void prepareCoupling(SuperLattice3D<T, DESCRIPTOR>& sLattice1,
                     SuperLattice3D<T, DESCRIPTOR>& sLattice2,
                     SuperGeometry3D<T>& superGeometry)
{
  OstreamManager clout( std::cout,"prepareCoupling" );
  clout << "Add lattice coupling" << endl;

  // (1) DEFINE THE REQUIRED COUPLINGS (POTENTIAL, FORCE, BOUNDARY)
  // Add the lattice couplings
  // The chemical potential coupling must come before the force coupling
  FreeEnergyChemicalPotentialGenerator3D<T, DESCRIPTOR> coupling1( alpha, kappa1, kappa2 );
  FreeEnergyForceGenerator3D<T, DESCRIPTOR> coupling2;
  // coupling for inlet BC, I.E. FOR MOVING WALLS
  FreeEnergyInletOutletGenerator3D<T, DESCRIPTOR> coupling3;

  // (2) ADD COUPLINGS TO LATTICE
  // The InletOutlet couplings must come after the Force coupling.
  sLattice1.addLatticeCoupling<DESCRIPTOR>( superGeometry, 1, coupling1, sLattice2 ); // INJECT POTENTIAL TO LATTICE 2
  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 1, coupling2, sLattice1 ); // INJECT FORCE TO LATTICE 1

  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 3, coupling3, sLattice1 ); // INJECT BOUNDARY TO LATTICE 1
  sLattice2.addLatticeCoupling<DESCRIPTOR>( superGeometry, 4, coupling3, sLattice1 ); // INJECT BOUNDARY TO LATTICE 1

  clout << "Add lattice coupling ... OK!" << endl;
}


//weitesgehend unberührt, bei vtk file mit ergänzter velocity Ausgabe
//evtl noch direkt jpeg Ausagbe des velocity Profils später einarbeiten
void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice1,
                 SuperLattice3D<T, DESCRIPTOR>& sLattice2,
                 int iT,
                 SuperGeometry3D<T>& superGeometry,
                 Timer<T>& timer,
                 UnitConverter<T, DESCRIPTOR> converter)
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "youngLaplace3d" );

  if ( iT==0 )
  {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice1, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice1 );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice1 );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }


  // Get statistics
  if ( iT%statIter==0 )
  {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice1.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice2.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 )
  {

    AnalyticalConst3D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperLatticeDensity3D<T, DESCRIPTOR> density1( sLattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity3D<T, DESCRIPTOR> density2( sLattice2 );
    density2.getName() = "phi";

    SuperIdentity3D<T,T> c1 (half*(density1+density2));
    c1.getName() = "density-fluid-1";
    SuperIdentity3D<T,T> c2 (half*(density1-density2));
    c2.getName() = "density-fluid-2";

    //neu für Implementierung der Geschwindigkeit in Paraview
    //je einzelne lattice Geschwindigkeiten + kombiniert
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity1( sLattice1, converter );
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity2( sLattice2, converter );
    SuperIdentity3D<T,T> velocityMittel (half*(velocity1+velocity2));
    velocityMittel.getName() = "velocity-Mittel";
    velocity1.getName() = "velocity1";
    velocity2.getName() = "velocity2";

    vtmWriter.addFunctor( density1 );
    vtmWriter.addFunctor( density2 );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );

    vtmWriter.addFunctor( velocity1 );
    vtmWriter.addFunctor( velocity2 );
    vtmWriter.addFunctor( velocityMittel );

    vtmWriter.write( iT );

    // calculate bulk pressure, pressure difference and surface tension
    if(iT%statIter==0)
    {
      AnalyticalConst3D<T,T> two_( 2. );
      AnalyticalConst3D<T,T> onefive_( 1.5 );
      AnalyticalConst3D<T,T> k1_( kappa1 );
      AnalyticalConst3D<T,T> k2_( kappa2 );
      AnalyticalConst3D<T,T> cs2_( 1./descriptors::invCs2<T,DESCRIPTOR>() );
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> two(two_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> onefive(onefive_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> k1(k1_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> k2(k2_, sLattice1);
      SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> cs2(cs2_, sLattice1);


      // DON'T NEED THAT RIGHT NOW

      // Calculation of bulk pressure:
      // c_1 = density of fluid 1; c_2 = density of fluid 2
      // p_bulk = rho*c_s^2 + kappa1 * (3/2*c_1^4 - 2*c_1^3 + 0.5*c_1^2)
      //                    + kappa2 * (3/2*c_2^4 - 2*c_2^3 + 0.5*c_2^2)
      // SuperIdentity3D<T,T> bulkPressure ( density1*cs2
      //         + k1*( onefive*c1*c1*c1*c1 - two*c1*c1*c1 + half*c1*c1 )
      //         + k2*( onefive*c2*c2*c2*c2 - two*c2*c2*c2 + half*c2*c2 ) );
      //
      // AnalyticalFfromSuperF3D<T, T> interpolPressure( bulkPressure, true, 1);
      // double position[3] = { 0.5*nx, 0.5*ny, 0.5*nz };
      // double pressureIn = 0.;
      // double pressureOut = 0.;
      // interpolPressure(&pressureIn, position);
      // position[0] = ((double)N/100.)*converter.getPhysDeltaX();
      // position[1] = ((double)N/100.)*converter.getPhysDeltaX();
      // position[2] = ((double)N/100.)*converter.getPhysDeltaX();
      // interpolPressure(&pressureOut, position);
      //
      // clout << "Pressure Difference: " << pressureIn-pressureOut << "  ;  ";
      // clout << "Surface Tension: " << radius*(pressureIn-pressureOut)/2 << std::endl;
      // clout << "Analytical Pressure Difference: " << alpha/(3.*radius) * (kappa1 + kappa2) << "  ;  ";
      clout << "Analytical Surface Tension: " << alpha/6. * (kappa1 + kappa2) << std::endl;
    }
  }
}


int main( int argc, char *argv[] )
{

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    (T)   N,      // resolution, default: N
    (T)   3.5,     // lattice relaxation time (tau)
    (T)   0.25*radius, //ny,     // charPhysLength: reference length of simulation geometry
    (T)   1./120.,//0.1,    // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__; default:1.e-6
    (T)   1.,      // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.      // physDensity: physical density in __kg / m^3__
  );
  // default Werte: N, 1., nx, 1.e-6, 0.1, 1.
  // H bzw. ny ist die charakteristische Länge (auch für Re)

  // Prints the converter log as console output
  converter.print();

  //Mobility Coefficient M für Peclet
  const T M = converter.getPhysDeltaT()*gama*(converter.getLatticeRelaxationTime()-0.5); //M=DeltaT*Gama*(Tau_g - 1/2)

  //AUSGABE wichtige Kennzahlen zur Überprüfung
  //physical dimensionsless parameter
  clout << "Reynolds Re=" << (2*(converter.getCharPhysVelocity())*radius*radius)/(ny*converter.getLatticeViscosity()) << std::endl; //Re=vL/nu=2v*a^2/H*nu
  clout << "Capillary Ca=" << (radius*2*(converter.getCharPhysVelocity())*converter.getLatticeViscosity())/((alpha/6. * (kappa1 + kappa2)) * ny) << std::endl; //Ca=mu*V/Sigma=a*2v*mu_c/H*Sigma
  clout << "Viscosity Ratio Lambda=1" << std::endl; //konst=1
  //numerical dimensionsless parameter
  clout << "Cahn Ch=" << alpha/radius << std::endl; //Ch=Sigma/r
  clout << "Peclet Pe*A=" << (2*(converter.getCharPhysVelocity())*radius*alpha)/(ny*M) << std::endl; //Pe=(2u*r*alpha)/(H*MA)

  // === 2nd Step: Prepare Geometry ===
  std::vector<T> extend = { nx, ny, nz };
  std::vector<T> origin = { 0, 0, 0 };
  IndicatorCuboid3D<T> cuboid( extend, origin );
#ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cGeometry( cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidGeometry3D<T> cGeometry( cuboid, converter.getPhysDeltaX() );
#endif

  // set periodic boundaries to the domain (x,y,z) -> hier nur x
  // nur x: true, false, false, heißt nur flow in x-Richtung
  cGeometry.setPeriodicity( true, false, false );

  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cGeometry );
  loadBalancer.print();

  // Instantiation of superGeometry
  SuperGeometry3D<T> superGeometry( cGeometry, loadBalancer );

  prepareGeometry( superGeometry, converter, cuboid);


  // === 3rd Step: Prepare Lattice ===
  SuperLattice3D<T, DESCRIPTOR> sLattice1( superGeometry );
  SuperLattice3D<T, DESCRIPTOR> sLattice2( superGeometry );

  ForcedBGKdynamics<T, DESCRIPTOR> bulkDynamics1 (
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,DESCRIPTOR>() );

  FreeEnergyBGKdynamics<T, DESCRIPTOR> bulkDynamics2 (
    converter.getLatticeRelaxationFrequency(), gama,
    instances::getBulkMomenta<T,DESCRIPTOR>() );

  //-------------------boundaries einbringen-----------------------------------
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sOnBC1( sLattice1 );
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sOnBC2( sLattice2 );
  createLocalBoundaryCondition3D<T, DESCRIPTOR> ( sOnBC1 );
  createLocalBoundaryCondition3D<T, DESCRIPTOR> ( sOnBC2 );

  prepareLattice( sLattice1, sLattice2, bulkDynamics1, bulkDynamics2, converter, sOnBC1, sOnBC2, superGeometry );

  prepareCoupling( sLattice1, sLattice2, superGeometry );

  SuperExternal3D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal1 ( superGeometry, sLattice1, sLattice1.getOverlap() );
  SuperExternal3D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal2 ( superGeometry, sLattice2, sLattice2.getOverlap() );

  // === 4th Step: Main Loop with Timer ===------------------------------------
  int iT = 0;
  clout << "starting simulation..." << endl;
  Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( iT=0; iT<=maxIter; ++iT ) {
    // Computation and output of the results
    getResults( sLattice1, sLattice2, iT, superGeometry, timer, converter );

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();

    // MPI communication for lattice data
    sLattice1.communicate();
    sLattice2.communicate();

    // Execute coupling between the two lattices
    sLattice1.executeCoupling();
    sExternal1.communicate();
    sExternal2.communicate();
    sLattice2.executeCoupling();
  }

  timer.stop();
  timer.printSummary();

}
