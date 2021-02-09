#pragma once

#include "algorithm.h"

#ifdef WINDOWS_VC
#define __DECLSPEC __declspec(dllexport)
#elif GNU_GCC
#define __DECLSPEC
#endif

//-----------------------------------------------------------------------------

// basic test for the interoperability of the library with the user level apps
// returns (numeric) '0xBABE'
extern "C" __DECLSPEC int sanityCheck();

// create/release a singleton for the pool of stiffness threads (within ins-
// tance threads) we are going to use. 'Release' destroys the pool object,
// i.e. you will need to 'Create' again, if you want to compute again. 'Clear'
// keeps the object but drops all instances (mesh threads). 'Create' returns
// 'SINGLE_POOL_ERROR' when trying to create more than 1 pool at a time and
// 'NO_ERROR' otherwise
extern "C" __DECLSPEC ErrorCode createForceDensityPool();
extern "C" __DECLSPEC void releaseForceDensityPool();
extern "C" __DECLSPEC void clearForceDensityPool();

// add a force-density instance to be run. The ID of the created instance is
// returned as second argument. The 'instance' represents a high level thread
// for a particular mesh. I.e. add as many instances as there are meshes you
// want to compute in parallel (unlimited, but keep reasonable, i.e. max 2-4).
// Returns 'BAD_INSTANCE_ERROR' for invalid pool objects (i.e. create the pool
// first) and 'NO_ERROR' otherwise. The number of 'stiffness' threads per mesh
// is internally limited to max 4 without error output. The 'ID' of the newly
// created instance is returned in the second argument
extern "C" __DECLSPEC ErrorCode addForceDensityInstance(
                               const int numStiffThreads,
                               int &instanceID);

// commit a point with specified coordinates/fixity to the specified instance.
// Returns 'BAD_INSTANCE_ERROR' for invalid pool objects, 'BAD_INDEX_ERROR'
// for instance ID's out of range and 'NO_ERROR' otherwise
extern "C" __DECLSPEC ErrorCode addInstancePointInfo(
                               const int instanceID,
                               const double x,
                               const double y,
                               const double z,
                               const bool fixed = false);

// alternative for supplying point loads different than default (i.e. 0.0).
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode addInstancePointInfoLoad(
                               const int instanceID,
                               const double x,
                               const double y,
                               const double z,
                               const double loadX,
                               const double loadY,
                               const double loadZ,
                               const bool fixed = false);

// commit a point with specified coordinates/fixity to the specified instance.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointInfo(
                               const int instanceID,
                               const int pointID,
                               const double x,
                               const double y,
                               const double z,
                               const bool fixed = false);

// convenience function for setting the X coordinate of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointX(
                               const int instanceID,
                               const int pointID,
                               const double x);

// convenience function for setting the Y coordinate of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointY(
                               const int instanceID,
                               const int pointID,
                               const double y);

// convenience function for setting the Z coordinate of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointZ(
                               const int instanceID,
                               const int pointID,
                               const double z);

// convenience function for setting the fixity value of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointFixed(
                               const int instanceID,
                               const int pointID,
                               const bool fixed);

// convenience function for setting the load triple for a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointLoad(
                               const int instanceID,
                               const int pointID,
                               const double loadX,
                               const double loadY,
                               const double loadZ);

// convenience function for setting the load value X for a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointLoadX(
                               const int instanceID,
                               const int pointID,
                               const double loadX);

// convenience function for setting the load value Y for a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointLoadY(
                               const int instanceID,
                               const int pointID,
                               const double loadY);

// convenience function for setting the load value Z for a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstancePointLoadZ(
                               const int instanceID,
                               const int pointID,
                               const double loadZ);

// read the input point data for the given instance and the given point index.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointInfo(
                               const int instanceID,
                               const int pointID,
                               double &x,
                               double &y,
                               double &z,
                               bool &fixed);

// convenience function for getting the X coordinate of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointX(
                               const int instanceID,
                               const int pointID,
                               double &x);

// convenience function for getting the Y coordinate of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointY(
                               const int instanceID,
                               const int pointID,
                               double &y);

// convenience function for getting the Z coordinate of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointZ(
                               const int instanceID,
                               const int pointID,
                               double &z);

// convenience function for getting the fixity value of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointFixed(
                               const int instanceID,
                               const int pointID,
                               bool &fixed);

// convenience function for getting the load values of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointLoad(
                               const int instanceID,
                               const int pointID,
                               double &loadX,
                               double &loadY,
                               double &loadZ);

// convenience function for getting the load X value of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointLoadX(
                               const int instanceID,
                               const int pointID,
                               double &loadX);

// convenience function for getting the load Y value of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointLoadY(
                               const int instanceID,
                               const int pointID,
                               double &loadY);

// convenience function for getting the load Z value of a particular point.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstancePointLoadZ(
                               const int instanceID,
                               const int pointID,
                               double &loadZ);

// commit an edge with specified endpoints/stiffness to the specified instance.
// Supply stiffness values for every thread, excess values (i.e. 4 stiffness
// values while only 2 stiffness threads created) are unused. The same error
// codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode addInstanceEdgeInfo(
                               const int instanceID,
                               const int ib,
                               const int ie,
                               const double stiffness1,
                               const double stiffness2,
                               const double stiffness3,
                               const double stiffness4);

// commit an edge with specified endpoints/stiffness to the specified instance.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstanceEdgeInfo(
                               const int instanceID,
                               const int edgeID,
                               const int ib,
                               const int ie,
                               const double stiff1,
                               const double stiff2,
                               const double stiff3,
                               const double stiff4);

// convenience function for setting the 'beginning' point of a particular edge.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstanceEdgeIB(
                               const int instanceID,
                               const int edgeID,
                               const int ib);

// convenience function for setting the 'ending' point of a particular edge.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstanceEdgeIE(
                               const int instanceID,
                               const int edgeIE,
                               const int ib);

// convenience function for setting one of the stiffness values of a particular
// edge. The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstanceEdgeStiffness(
                               const int instanceID,
                               const int edgeID,
                               const int threadID,
                               const double stiffness);

// convenience function for setting all 4 stiffness values at once. I.e. one
// for each thread (max. 4 stiffness threads per instance). The same error
// codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstanceEdgeStiffnessQuad(
                               const int instanceID,
                               const int edgeID,
                               const double stiffness1,
                               const double stiffness2,
                               const double stiffness3,
                               const double stiffness4);

// read the input edge data for the given instance and the given edge index.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceEdgeInfo(
                               const int instanceID,
                               const int edgeID,
                               int &ib,
                               int &ie,
                               double &stiff1,
                               double &stiff2,
                               double &stiff3,
                               double &stiff4);

// convenience function for getting the 'beginning' point of a particular edge.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceEdgeIB(
                               const int instanceID,
                               const int edgeID,
                               int &ib);

// convenience function for getting the 'ending' point of a particular edge.
// The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceEdgeIE(
                               const int instanceID,
                               const int edgeID,
                               int &ie);

// convenience function for getting one of the stiffness values of a particular
// edge. The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceEdgeStiffness(
                               const int instanceID,
                               const int edgeID,
                               const int threadID,
                               double &stiffness);

// convenience function for getting all 4 stiffness values at once. I.e. one
// for each thread (max. 4 stiffness threads per instance). The same error
// codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceEdgeStiffnessQuad(
                               const int instanceID,
                               const int edgeID,
                               double &stiffness1,
                               double &stiffness2,
                               double &stiffness3,
                               double &stiffness4);

// allow assigning 'uniform' (i.e. constant) load to each thread. This allows
// a multi-load computation for the same (mesh) instance. The same error codes
// as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode setInstanceUniformLoad(
                               const int instanceID,
                               const int threadID,
                               const double loadX,
                               const double loadY,
                               const double loadZ);

// allow reading 'uniform' (i.e. constant) loads for each thread. This allows
// a multi-load computation for the same (mesh) instance. The same error codes
// as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceUniformLoad(
                               const int instanceID,
                               const int threadID,
                               double &loadX,
                               double &loadY,
                               double &loadZ);

// allow reading result data for the specified mesh ('instanceID'), stiffness-
// thread ('threadID'), and specified point ('pointID'). Coordinates only,
// loads are not supposed to change, i.e. remain the same as the input. The
// same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceResultPoint(
                               const int instanceID,
                               const int threadID,
                               const int pointID,
                               double &x,
                               double &y,
                               double &z);

// convenience function for getting the X coordinate of the specified point
// for the specified instance/thread. Error codes see 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceResultPointX(
                               const int instanceID,
                               const int threadID,
                               const int pointID,
                               double &x);

// convenience function for getting the Y coordinate of the specified point
// for the specified instance/thread. Error codes see 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceResultPointY(
                               const int instanceID,
                               const int threadID,
                               const int pointID,
                               double &y);

// convenience function for getting the Z coordinate of the specified point
// for the specified instance/thread. Error codes see 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode getInstanceResultPointZ(
                               const int instanceID,
                               const int threadID,
                               const int pointID,
                               double &z);

// return the number of points/edges/instances. Returns -1 in case of errors
// (ID out of bounds, released instance, etc) or a valid count if there are
// no errors
extern "C" __DECLSPEC int getNumForceDensityInstances();
extern "C" __DECLSPEC int getNumInstancePoints(const int instanceID);
extern "C" __DECLSPEC int getNumInstanceEdges(const int instanceID);
extern "C" __DECLSPEC int getNumInstanceThreads(const int instanceID);

// enable/disable the entire force-density instance (i.e. top level threads
// for particular meshes). This might be convenient for iterative processing,
// with one mesh requiring more/less iterations than the other. Enabled by
// default. The same error codes as 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode enableForceDensityInstance(
                               const int instanceID,
                               const bool enable);

// allow enabling/disabling separate stiffness threads for a given instance.
// This might also be convenient, all threads are enabled by default. For error
// codes see 'addInstancePointInfo'
extern "C" __DECLSPEC ErrorCode enableInstanceThread(
                               const int instanceID,
                               const int threadID,
                               const bool enable);

// run all collected/added force-density instances (mesh computation) at once
// in separate threads, which in turn initiate further threads (multiple stiff-
// ness computations per mesh)
extern "C" __DECLSPEC ErrorCode runForceDensityInstances();


// OLD STUFF, dont use!

// additional parameters for the calculator. Allow the surface to be influen-
// ced by external loads. Interpret as a non-normalized 3d vector with compo-
// nents representing loads in a particular direction
extern "C" __DECLSPEC ErrorCode setUniformSurfaceLoad(const double,
                                                      const double,
                                                      const double);

extern "C" __DECLSPEC ErrorCode getUniformSurfaceLoad(double&,
                                                      double&,
                                                      double&);

