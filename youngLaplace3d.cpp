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
const int N  = 100;
const T nx   = 100.;
const T radius = 0.25 * nx;
const T alpha = 1.5;     // Interfacial width         [lattice units]
const T kappa1 = 0.0075; // For surface tensions      [lattice units]
const T kappa2 = 0.005;  // For surface tensions      [lattice units]
const T gama = 1.;       // For mobility of interface [lattice units]

//h_i neu (analog microFluidics2d) zur Einführung von freeEnergyBoundary
//h_i= Parameter related to resulting contact angle of the boundary. [lattice units]
const T h1 = 0.;                  // Contact angle 90 degrees   [lattice units]
const T h2 = 0.;                  // Contact angle 90 degrees   [lattice units]


//neu eingeführt für Geometry MN, dabei entsprechen sie den Raumkoord. nx
const T lengthX = 100.;
const T lengthY = 100.;
const T lengthZ = 100.;

const int maxIter  = 1000; //default 60.000 -> hier Simulationsschritte einstellen
const int vtkIter  = 200;
const int statIter = 200;


void prepareGeometry( SuperGeometry3D<T>& superGeometry,
                      UnitConverter<T,DESCRIPTOR> const& converter,
                      IndicatorF3D<T>& indicator)
{

  T eps = converter.getConversionFactorLength();

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  //default, alles 1: superGeometry.rename( 0,1 );

  //---------------------MATERIAL NUMBERS --------------------------------------
  Vector<T,3> extendGeometryInOut( lengthX,lengthY,lengthZ); //erstmal alle Längen = nx
  Vector<T,3> origin;
  IndicatorCuboid3D<T> cuboid(extendGeometryInOut, origin);


  //replace MN=0 mit 2 (2=Wand) für cuboid
  //NEU: mit Cuboid Indicator von (0,0,0) auf (100,100,100)
  superGeometry.rename(0, 2);

  //ersetze 2(boundary) mit 1(fluid) mit offset 1; ändere MN von 2 auf 1 für
  //innere Zellen mit der Form "void rename (fromM, toN, offsetX, offsetY, offsetZ)" da 3D
  superGeometry.rename(2, 1, 1, 1, 1);
  /*normalerweise extend=(nx,nx,nx), hier will ich aber über (0,nx,nx) die Fläche
  für Inflow und Outflow aufspannen, von origin(0,0,0) und origin(nx,0,0), daher
  für jede MN einzeln neu definieren. Da 0 nicht als Wert für Vektor erlaubt,
  wird mit eps gearbeitet. */

  // Ändere MN=3 Inflow
  Vector<T,3> origin3 (-eps, -eps, -eps);
  Vector<T,3> extendGeometryInOut3( +2*eps, lengthY+2*eps, lengthZ+2*eps);
  IndicatorCuboid3D<T> inflow(extendGeometryInOut3, origin3);
  //superGeometry.rename(2, 3, 1,inflow);
  superGeometry.rename(1, 3, inflow);
  superGeometry.rename(2, 3, inflow); //void rename (from, to, fluidMN, indicator functor condition)

  //Ändere MN=4 Outflow
  Vector<T,3> origin4(lengthX-eps, -eps, -eps);
  Vector<T,3> extendGeometryInOut4( lengthX+2*eps, lengthY+2*eps, lengthZ+2*eps);
  IndicatorCuboid3D<T> outflow(extendGeometryInOut4, origin4);
  //superGeometry.rename(2, 4, 1, outflow);
  superGeometry.rename(1, 4, outflow);
  superGeometry.rename(2, 4, outflow);
  // numeric_limits<T>::epsilon)()

  //obere Wand MN=5 mit Wandgeschwindigkeit
  Vector<T,3> origin5(-eps, 100.-eps, -eps);
  Vector<T,3> extendGeometryInOut5( lengthX+2*eps, +2*eps, lengthZ+2*eps);
  IndicatorCuboid3D<T> oben(extendGeometryInOut5, origin5);
  superGeometry.rename(2, 5, 1, oben);

  //noch mit unterer Wand mit Wandgeschwindigkeit analog MN=6!
  Vector<T,3> origin6(-eps, -eps, -eps);
  Vector<T,3> extendGeometryInOut6( lengthX+2*eps, +2*eps, lengthZ+2*eps);
  IndicatorCuboid3D<T> unten(extendGeometryInOut6, origin6);
  superGeometry.rename(2, 6, 1, unten);

  //----------------------------------------------------------------------------
  // Removes all not needed boundary voxels outside the surface
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
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBCvel,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBCvel2,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBCvel3,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBCvel4,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBCenergy1,
                     sOnLatticeBoundaryCondition3D<T,DESCRIPTOR>& sOnBCenergy2,
                     SuperGeometry3D<T>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  //neu für omega
  const T omega = converter.getLatticeRelaxationFrequency();

  //define lattice Dynamics, jeweils beide Lattices beachten
  //bulkDynamics1 = forced ForcedBGKdynamics
  //bulkDynamics2 = FreeEnergyBGKdynamics

  //MN=0 -> no dynamics
  sLattice1.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 0, &instances::getNoDynamics<T, DESCRIPTOR>() );

  //MN=1 Fluid hat bulkdynamics
  sLattice1.defineDynamics( superGeometry, 1, &bulkDynamics1 );
  sLattice2.defineDynamics( superGeometry, 1, &bulkDynamics2 );

  //MN=2
  //geplant: Wände erstmal bounce back, jetzt: NotDynamics -> über BC FreeEnergyWallBoundary
  sLattice1.defineDynamics( superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );
  sLattice2.defineDynamics( superGeometry, 2, &instances::getNoDynamics<T, DESCRIPTOR>() );

  //MN=3,4 in/outflow erstmal bulk dynamics
  sLattice1.defineDynamics(superGeometry, 3, &bulkDynamics1);
  sLattice2.defineDynamics(superGeometry, 3, &bulkDynamics2);

  sLattice1.defineDynamics(superGeometry, 4, &bulkDynamics1);
  sLattice2.defineDynamics(superGeometry, 4, &bulkDynamics2);

  //MN=5 obere Wand, eigene Dynamics wegen Wandgeschwindigkeit
  sLattice1.defineDynamics(superGeometry, 5, &bulkDynamics1);
  sLattice2.defineDynamics(superGeometry, 5, &bulkDynamics2);

  // MN=6 untere Wandgeschwindigkeit
  sLattice1.defineDynamics(superGeometry, 6, &bulkDynamics1);
  sLattice2.defineDynamics(superGeometry, 6, &bulkDynamics2);

  //--------------------------BOUNDARIES----------------------------------------
  //add velocity boundary für bewegende Wand MN=5
  sOnBCvel.addVelocityBoundary( superGeometry, 5, omega );
  sOnBCvel2.addVelocityBoundary( superGeometry, 5, omega );
  sOnBCvel3.addVelocityBoundary( superGeometry, 6, omega );
  sOnBCvel4.addVelocityBoundary( superGeometry, 6, omega );
  //add free Energy boundaries für Wände MN=2
  //gemäß Doxygen für 3D, 2Phasen: x(superGeometry, MN, alpha, kappa1, kappa2, h1, h2, latticeNumber)
  sOnBCenergy1.addFreeEnergyWallBoundary( superGeometry, 2, alpha, kappa1, kappa2, h1, h2, 1 );
  sOnBCenergy2.addFreeEnergyWallBoundary( superGeometry, 2, alpha, kappa1, kappa2, h1, h2, 2 );

  //----------------------------------------------------------------------------
  AnalyticalConst3D<T,T> rhoF( T( 1 ) ); //setze rho=1
  AnalyticalConst3D<T,T> uF( T( 0 ), T( 0 ), T( 0 ) ); //setze vec(v)=0=(0,0,0)

  //hier aus default Code prinipiell selbes Spiel wie oben mit rho=1, v=0. Aber mit phi für
  //bulk initial conditions
  //rho=1 und v=0 hier und oben quasi doppelt definiert -> kann später vereinfacht werden

  std::vector<T> v( 3,T() );
  AnalyticalConst3D<T,T> zeroVelocity( v ); // Geschwindigkeit (eigentlich schon oben festgelegt)
  AnalyticalConst3D<T,T> zero ( 0. ); // Null
  AnalyticalConst3D<T,T> one ( 1. ); // Eins
  SmoothIndicatorSphere3D<T,T> sphere( {nx/2., nx/2., nx/2.}, radius, 10.*alpha ); //Tropfen
  AnalyticalIdentity3D<T,T> rho( one ); //rho=1
  AnalyticalIdentity3D<T,T> phi( one - sphere - sphere );

  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 2, 3, 4, 5, 6});

  //bei sLattice 1 mit rho, sLattice2 mit phi
  //(bulkIndicator, rhoF, uF) ersetzt (superGeometry, 1, rho, zeroVelocity)
  sLattice1.iniEquilibrium( bulkIndicator, rhoF, uF );
  sLattice2.iniEquilibrium( bulkIndicator, phi, uF );

  sLattice1.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice2.defineRhoU( bulkIndicator, phi, uF );

  //anstatt über converter direkt angeben
  clout << converter.getCharLatticeVelocity() << std::endl;
  AnalyticalConst3D<T,T> uTop( converter.getCharLatticeVelocity(), T( 0 ), T( 0 ) );
  AnalyticalConst3D<T,T> uDown( -converter.getCharLatticeVelocity(), T( 0 ), T( 0 ) );
  //Alternative Überlegung ohne converter. Weil: getCharLatticeVelocity holt Wandgeschwindigkeit
  //aus charPhysVelocity aus den Werten der main Funktion. Die entsprechen jedoch der
  //Strömungsgeschwnidigkeit. So wäre also Wandgeschwindigkeit=Strömungsgeschwnidigkeit?
  //hier Geschwindigkeit anstatt über converter, direkt angeben, Annahme: uTop(v_x, v_y, v_z)
  //v_x=1000 m/s entspräche:
  //clout << T(1000.) << std::endl;
  //AnalyticalConst3D<T,T> uTop( T(1000.) , T( 0 ), T( 0 ) );

  //MN=5 ist obere Wandgeschw. uTop
  sLattice1.defineU( superGeometry, 5, uTop );
  sLattice2.defineU( superGeometry, 5, uTop );
  sLattice1.defineU( superGeometry, 6, uDown );
  sLattice2.defineU( superGeometry, 6, uDown );

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


//Kopplung erstmal unberührt lassen
void prepareCoupling(SuperLattice3D<T, DESCRIPTOR>& sLattice1,
                     SuperLattice3D<T, DESCRIPTOR>& sLattice2)
{

  OstreamManager clout( std::cout,"prepareCoupling" );
  clout << "Add lattice coupling" << endl;

  // Add the lattice couplings
  // The chemical potential coupling must come before the force coupling
  FreeEnergyChemicalPotentialGenerator3D<T, DESCRIPTOR> coupling1(
    alpha, kappa1, kappa2);
  FreeEnergyForceGenerator3D<T, DESCRIPTOR> coupling2;

  sLattice1.addLatticeCoupling( coupling1, sLattice2 );
  sLattice2.addLatticeCoupling( coupling2, sLattice1 );

  clout << "Add lattice coupling ... OK!" << endl;
}


//weitesgehend unberührt, bei vtk file mit velocity Ausgabe
//evtl noch direkt jpeg Ausagbe des velocity Profils später einarbeiten
void getResults( SuperLattice3D<T, DESCRIPTOR>& sLattice2,
                 SuperLattice3D<T, DESCRIPTOR>& sLattice1,
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

      // Calculation of bulk pressure:
      // c_1 = density of fluid 1; c_2 = density of fluid 2
      // p_bulk = rho*c_s^2 + kappa1 * (3/2*c_1^4 - 2*c_1^3 + 0.5*c_1^2)
      //                    + kappa2 * (3/2*c_2^4 - 2*c_2^3 + 0.5*c_2^2)
      SuperIdentity3D<T,T> bulkPressure ( density1*cs2
              + k1*( onefive*c1*c1*c1*c1 - two*c1*c1*c1 + half*c1*c1 )
              + k2*( onefive*c2*c2*c2*c2 - two*c2*c2*c2 + half*c2*c2 ) );

      AnalyticalFfromSuperF3D<T, T> interpolPressure( bulkPressure, true, 1);
      double position[3] = { 0.5*nx, 0.5*nx, 0.5*nx };
      double pressureIn = 0.;
      double pressureOut = 0.;
      interpolPressure(&pressureIn, position);
      position[0] = ((double)N/100.)*converter.getPhysDeltaX();
      position[1] = ((double)N/100.)*converter.getPhysDeltaX();
      position[2] = ((double)N/100.)*converter.getPhysDeltaX();
      interpolPressure(&pressureOut, position);

      clout << "Pressure Difference: " << pressureIn-pressureOut << "  ;  ";
      clout << "Surface Tension: " << radius*(pressureIn-pressureOut)/2 << std::endl;
      clout << "Analytical Pressure Difference: " << alpha/(3.*radius) * (kappa1 + kappa2) << "  ;  ";
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
    (T)   N,      // resolution
    (T)   1.,     // lattice relaxation time (tau)
    (T)   nx,     // charPhysLength: reference length of simulation geometry
    (T)   0.1,  // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__; default:1.e-6
    (T)   0.1,    // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.      // physDensity: physical density in __kg / m^3__
  );
  //default Werte: N, 1., nx, 1.e-6, 0.1, 1.

  // Prints the converter log as console output
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  std::vector<T> extend = { nx, nx, nx };
  std::vector<T> origin = { 0, 0, 0 };
  IndicatorCuboid3D<T> cuboid(extend,origin);
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
  SuperGeometry3D<T> superGeometry( cGeometry,loadBalancer );

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

  //velocity boundary
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sOnBCvel(sLattice1);
  createLocalBoundaryCondition3D<T,DESCRIPTOR>( sOnBCvel );
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sOnBCvel2(sLattice2);
  createLocalBoundaryCondition3D<T,DESCRIPTOR>( sOnBCvel2 );
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sOnBCvel3(sLattice1);
  createLocalBoundaryCondition3D<T,DESCRIPTOR>( sOnBCvel3 );
  sOnLatticeBoundaryCondition3D<T,DESCRIPTOR> sOnBCvel4(sLattice1);
  createLocalBoundaryCondition3D<T,DESCRIPTOR>( sOnBCvel4 );
  //davor als createInterpBoundaryCondition3D, jetzt local
  //createInterpBoundaryCondition3D<T,DESCRIPTOR, ForcedBGKdynamics<T, DESCRIPTOR> >( sOnBCvel );

  //freeEnergy boundaries
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sOnBCenergy1( sLattice1 );
  sOnLatticeBoundaryCondition3D<T, DESCRIPTOR> sOnBCenergy2( sLattice2 );
  createLocalBoundaryCondition3D<T, DESCRIPTOR> (sOnBCenergy1);
  createLocalBoundaryCondition3D<T, DESCRIPTOR> (sOnBCenergy2);

  prepareLattice(sLattice1, sLattice2, bulkDynamics1, bulkDynamics2, converter, sOnBCvel, sOnBCvel2, sOnBCvel3, sOnBCvel4, sOnBCenergy1, sOnBCenergy2, superGeometry );

  prepareCoupling( sLattice1, sLattice2);

  SuperExternal3D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal1 (superGeometry, sLattice1, sLattice1.getOverlap() );
  SuperExternal3D<T,DESCRIPTOR,CHEM_POTENTIAL> sExternal2 (superGeometry, sLattice2, sLattice2.getOverlap() );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << endl;
  Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( iT=0; iT<=maxIter; ++iT ) {
    // Computation and output of the results
    getResults( sLattice2, sLattice1, iT, superGeometry, timer, converter );

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
