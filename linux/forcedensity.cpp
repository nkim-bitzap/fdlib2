#include "forcedensity.h"
#include "exception.h"

ForceDensityPool *pool = 0;

//-----------------------------------------------------------------------------

int sanityCheck() {
  return 0xBABE;
}

//-----------------------------------------------------------------------------

ErrorCode createForceDensityPool() {
  if (pool) return SINGLE_POOL_ERROR;
  else {
    pool = new ForceDensityPool();
    return NO_ERROR;
  }
}

//-----------------------------------------------------------------------------

void clearForceDensityPool() {
  if (pool) pool->releaseInstances();
}

//-----------------------------------------------------------------------------

void releaseForceDensityPool() {
  if (pool) {
    pool->releaseInstances();

    delete pool;
    pool = 0;
  }
}

//-----------------------------------------------------------------------------

ErrorCode addForceDensityInstance(const int numThreads, int &instanceID) {
  if (!pool) return BAD_INSTANCE_ERROR;

  instanceID = pool->getNumInstances();
  pool->addInstance(new ForceDensity(numThreads));
  return NO_ERROR;
}

//-----------------------------------------------------------------------------

int getNumForceDensityInstances() {
  if (pool) return pool->getNumInstances();
  else return -1;
}

//-----------------------------------------------------------------------------

ErrorCode addInstancePointInfo(const int instanceID,
                               const double x,
                               const double y,
                               const double z,
                               const bool fixed)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->addPoint(x, y, z, fixed);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode addInstancePointInfoLoad(const int instanceID,
                                   const double x,
                                   const double y,
                                   const double z,
                                   const double loadX,
                                   const double loadY,
                                   const double loadZ,
                                   const bool fixed)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->addPoint(
      x, y, z, LoadInfo(loadX, loadY, loadZ), fixed);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointInfo(const int instanceID,
                               const int pointID,
                               const double x,
                               const double y,
                               const double z,
                               const bool fixed)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->setPoint(
      pointID, PointInfo(x, y, z, fixed));
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointX(const int instanceID,
                            const int pointID,
                            const double x)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getPoint(pointID).setX(x);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointY(const int instanceID,
                            const int pointID,
                            const double y)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getPoint(pointID).setY(y);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointZ(const int instanceID,
                            const int pointID,
                            const double z)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getPoint(pointID).setZ(z);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointFixed(const int instanceID,
                            const int pointID,
                            const bool fixed)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getPoint(pointID).setFixed(fixed);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointLoad(const int instanceID,
                               const int pointID,
                               const double loadX,
                               const double loadY,
                               const double loadZ)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getPoint(
      pointID).setLoad(loadX, loadY, loadZ);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointLoadX(const int instanceID,
                                const int pointID,
                                const double loadX)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getPoint(
      pointID).getLoadInfo().setLoadX(loadX);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointLoadY(const int instanceID,
                                const int pointID,
                                const double loadY)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getPoint(
      pointID).getLoadInfo().setLoadY(loadY);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstancePointLoadZ(const int instanceID,
                                const int pointID,
                                const double loadZ)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getPoint(
      pointID).getLoadInfo().setLoadZ(loadZ);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointInfo(const int instanceID,
                               const int pointID,
                               double &x,
                               double &y,
                               double &z,
                               bool &fixed)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    const PointInfo &pi =
      pool->getInstance(instanceID)->getPoint(pointID);

    x = pi.x();
    y = pi.y();
    z = pi.z();

    fixed = pi.getFixed();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointX(const int instanceID,
                            const int pointID,
                            double &x)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    x = pool->getInstance(instanceID)->getPoint(pointID).x();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointY(const int instanceID,
                            const int pointID,
                            double &y)
{
  if (!pool) return BAD_INSTANCE_ERROR;
  
  try {
    y = pool->getInstance(instanceID)->getPoint(pointID).y();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }
  
  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointZ(const int instanceID,
                            const int pointID,
                            double &z)
{
  if (!pool) return BAD_INSTANCE_ERROR;
  
  try {
    z = pool->getInstance(instanceID)->getPoint(pointID).z();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }
  
  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointFixed(const int instanceID,
                                const int pointID,
                                bool &fixed)
{
  if (!pool) return BAD_INSTANCE_ERROR;
  
  try {
    fixed =
      pool->getInstance(instanceID)->getPoint(pointID).getFixed();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }
  
  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointLoad(const int instanceID,
                               const int pointID,
                               double &loadX,
                               double &loadY,
                               double &loadZ)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    const LoadInfo &li = pool->getInstance(instanceID)->
      getPoint(pointID).getLoadInfo();

    loadX = li.loadX();
    loadY = li.loadY();
    loadZ = li.loadZ();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointLoadX(const int instanceID,
                                const int pointID,
                                double &loadX)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    const LoadInfo &li = pool->getInstance(instanceID)->
      getPoint(pointID).getLoadInfo();

    loadX = li.loadX();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointLoadY(const int instanceID,
                                const int pointID,
                                double &loadY)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    const LoadInfo &li = pool->getInstance(instanceID)->
      getPoint(pointID).getLoadInfo();

    loadY = li.loadY();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstancePointLoadZ(const int instanceID,
                                const int pointID,
                                double &loadZ)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    const LoadInfo &li = pool->getInstance(instanceID)->
      getPoint(pointID).getLoadInfo();

    loadZ = li.loadZ();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode addInstanceEdgeInfo(const int instanceID,
                              const int ib,
                              const int ie,
                              const double stiff1,
                              const double stiff2,
                              const double stiff3,
                              const double stiff4)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->addEdge(
      ib, ie, StiffnessQuad(stiff1, stiff2, stiff3, stiff4));
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstanceEdgeInfo(const int instanceID,
                              const int edgeID,
                              const int ib,
                              const int ie,
                              const double stiff1,
                              const double stiff2,
                              const double stiff3,
                              const double stiff4)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->setEdge(
      edgeID, EdgeInfo(
        ib, ie, StiffnessQuad(stiff1, stiff2, stiff3, stiff4)));
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstanceEdgeIB(const int instanceID,
                            const int edgeID,
                            const int ib)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getEdge(edgeID).setIB(ib);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstanceEdgeIE(const int instanceID,
                            const int edgeID,
                            const int ie)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getEdge(edgeID).setIE(ie);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstanceEdgeStiffness(const int instanceID,
                                   const int edgeID,
                                   const int threadID,
                                   const double stiffness)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->
     getEdge(edgeID).getStiffnessQuad().setStiffness(
       threadID, stiffness);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstanceEdgeStiffnessQuad(const int instanceID,
                                       const int edgeID,
                                       const double stiffness1,
                                       const double stiffness2,
                                       const double stiffness3,
                                       const double stiffness4)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->getEdge(edgeID).setStiffness(
      stiffness1, stiffness2, stiffness3, stiffness4);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode setInstanceUniformLoad(const int instanceID,
                                 const int threadID,
                                 const double loadX,
                                 const double loadY,
                                 const double loadZ)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->setUniformLoad(
      threadID, loadX, loadY, loadZ);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceUniformLoad(const int instanceID,
                                 const int threadID,
                                 double &loadX,
                                 double &loadY,
                                 double &loadZ)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    const LoadInfo &li =
      pool->getInstance(instanceID)->getUniformLoad(threadID);

    loadX = li.loadX();
    loadY = li.loadY();
    loadZ = li.loadZ();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceEdgeInfo(const int instanceID,
                              const int edgeID,
                              int &ib,
                              int &ie,
                              double &stiff1,
                              double &stiff2,
                              double &stiff3,
                              double &stiff4)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    const EdgeInfo &ei =
      pool->getInstance(instanceID)->getEdge(edgeID);

    ib = ei.ib();
    ie = ei.ie();

    stiff1 = ei.getStiffnessQuad().getStiffness(0);
    stiff2 = ei.getStiffnessQuad().getStiffness(1);
    stiff3 = ei.getStiffnessQuad().getStiffness(2);
    stiff4 = ei.getStiffnessQuad().getStiffness(3);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceEdgeIB(const int instanceID,
                            const int edgeID,
                            int &ib)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    ib =
      pool->getInstance(instanceID)->getEdge(edgeID).ib();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceEdgeIE(const int instanceID,
                            const int edgeID,
                            int &ie)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    ie =
      pool->getInstance(instanceID)->getEdge(edgeID).ie();
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceEdgeStiffness(const int instanceID,
                                   const int edgeID,
                                   const int threadID,
                                   double &stiffness)
{
  if (!pool) return BAD_INSTANCE_ERROR;
  
  try {
    stiffness = pool->getInstance(instanceID)->
      getEdge(edgeID).getStiffnessQuad().getStiffness(threadID);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceEdgeStiffnessQuad(const int instanceID,
                                       const int edgeID,
                                       double &stiffness1,
                                       double &stiffness2,
                                       double &stiffness3,
                                       double &stiffness4)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    const StiffnessQuad &sq =pool->getInstance(instanceID)->
      getEdge(edgeID).getStiffnessQuad();

    stiffness1 = sq.getStiffness(0);
    stiffness2 = sq.getStiffness(1);
    stiffness3 = sq.getStiffness(2);
    stiffness4 = sq.getStiffness(3);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }
  
  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceResultPoint(const int instanceID,
                                 const int threadID,
                                 const int pointID,
                                 double &x,
                                 double &y,
                                 double &z)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    ForceDensity *fd = pool->getInstance(instanceID);
    const Eigen::MatrixXd &result = fd->getResult(threadID);

    x = result(pointID, 0);
    y = result(pointID, 1);
    z = result(pointID, 2);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceResultPointX(const int instanceID,
                                  const int threadID,
                                  const int pointID,
                                  double &x)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    ForceDensity *fd = pool->getInstance(instanceID);
    const Eigen::MatrixXd &result = fd->getResult(threadID);

    x = result(pointID, 0);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceResultPointY(const int instanceID,
                                  const int threadID,
                                  const int pointID,
                                  double &y)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    ForceDensity *fd = pool->getInstance(instanceID);
    const Eigen::MatrixXd &result = fd->getResult(threadID);

    y = result(pointID, 1);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode getInstanceResultPointZ(const int instanceID,
                                  const int threadID,
                                  const int pointID,
                                  double &z)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    ForceDensity *fd = pool->getInstance(instanceID);
    const Eigen::MatrixXd &result = fd->getResult(threadID);

    z = result(pointID, 2);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

int getNumInstancePoints(const int instanceID) {
  if (!pool) return -1;

  try {
    return pool->getInstance(instanceID)->getPoints().size();
  }
  catch (IndexAccessException) {
    return -1;
  }
}

//-----------------------------------------------------------------------------

int getNumInstanceEdges(const int instanceID) {
  if (!pool) return -1;

  try {
    return pool->getInstance(instanceID)->getEdges().size();
  }
  catch (IndexAccessException) {
    return -1;
  }
}

//-----------------------------------------------------------------------------

int getNumInstanceThreads(const int instanceID) {
  if (!pool) return -1;

  try {
    return pool->getInstance(instanceID)->getNumStiffnessThreads();
  }
  catch (IndexAccessException) {
    return -1;
  }
}

//-----------------------------------------------------------------------------

ErrorCode enableForceDensityInstance(const int instanceID,
                                     const bool enable)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->enableInstance(enable);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode enableInstanceThread(const int instanceID,
                               const int threadID,
                               const bool enable)
{
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->getInstance(instanceID)->
      enableStiffnessThread(threadID, enable);
  }
  catch (IndexAccessException) {
    return BAD_INDEX_ERROR;
  }

  return NO_ERROR;
}

//-----------------------------------------------------------------------------

ErrorCode runForceDensityInstances() {
  if (!pool) return BAD_INSTANCE_ERROR;

  try {
    pool->runInstances();
  }
  catch (const BadInputException &e) {
    return EMPTY_ALGORITHM_INPUT_ERROR;
  }
  catch (const BadEdgeException &e) {
    return MALFORMED_EDGE_ERROR;
  }
  catch (const IndexAccessException &e) {
    return BAD_INDEX_ERROR;
  }
  catch (const SolverComputeException &e) {
    return BAD_EIGEN_COMPUTE_ERROR;
  }
  catch (const SolverException &e) {
    return BAD_EIGEN_SOLVE_ERROR;
  }
  catch (...) {
    return UNKNOWN_ERROR;
  }

  return NO_ERROR;
}

/*
//-----------------------------------------------------------------------------

ErrorCode setUniformSurfaceLoad(const double x,
                                const double y,
                                const double z)
{
  if (!_fd) return BAD_INSTANCE_ERROR;
  else {
    _fd->setUniformSurfaceLoad(x, y, z);
    return NO_ERROR;
  }
}

//-----------------------------------------------------------------------------

ErrorCode getUniformSurfaceLoad(double &x,
                                double &y,
                                double &z)
{
  if (!_fd) return BAD_INSTANCE_ERROR;
  else {
    _fd->getUniformSurfaceLoad(x, y, z);
    return NO_ERROR;
  }
}
*/
