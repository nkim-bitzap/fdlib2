#include "algorithm.h"
#include "exception.h"

#include <iostream>
#include <list>
#include <thread>

//-----------------------------------------------------------------------------
// LoadInfo implementation
//-----------------------------------------------------------------------------

LoadInfo::LoadInfo()
: _loadX(0.0),
  _loadY(0.0),
  _loadZ(0.0)
{}

//-----------------------------------------------------------------------------

LoadInfo::LoadInfo(double x, double y, double z)
: _loadX(x),
  _loadY(y),
  _loadZ(z)
{}

//-----------------------------------------------------------------------------

double LoadInfo::loadX() const {
  return _loadX;
}

//-----------------------------------------------------------------------------

double LoadInfo::loadY() const {
  return _loadY;
}

//-----------------------------------------------------------------------------

double LoadInfo::loadZ() const {
  return _loadZ;
}

//-----------------------------------------------------------------------------

void LoadInfo::setLoadX(double x) {
  _loadX = x;
}

//-----------------------------------------------------------------------------

void LoadInfo::setLoadY(double y) {
  _loadY = y;
}

//-----------------------------------------------------------------------------

void LoadInfo::setLoadZ(double z) {
  _loadZ = z;
}
    
//-----------------------------------------------------------------------------

void LoadInfo::getLoad(double &x, double &y, double &z) {
  x = _loadX;
  y = _loadY;
  z = _loadZ;
}

//-----------------------------------------------------------------------------

void LoadInfo::setLoad(double x, double y, double z) {
  _loadX = x;
  _loadY = y;
  _loadZ = z;
}

//-----------------------------------------------------------------------------
// StiffnessQuad implementation
//-----------------------------------------------------------------------------

StiffnessQuad::StiffnessQuad()
: _stiffness(std::vector<double>(MAX_NUM_STIFFNESS_THREADS, 0.0))
{}

//-----------------------------------------------------------------------------

StiffnessQuad::StiffnessQuad(double stiffness)
: _stiffness(std::vector<double>(MAX_NUM_STIFFNESS_THREADS, 0.0))
{
  _stiffness[0] = stiffness;
}

//-----------------------------------------------------------------------------
    
StiffnessQuad::StiffnessQuad(double stiffness1,
                             double stiffness2)
: _stiffness(std::vector<double>(MAX_NUM_STIFFNESS_THREADS, 0.0))
{
  _stiffness[0] = stiffness1;
  _stiffness[1] = stiffness2;
}
    
//-----------------------------------------------------------------------------

StiffnessQuad::StiffnessQuad(double stiffness1,
                             double stiffness2,
                             double stiffness3)
: _stiffness(std::vector<double>(MAX_NUM_STIFFNESS_THREADS, 0.0))
{
  _stiffness[0] = stiffness1;
  _stiffness[1] = stiffness2;
  _stiffness[2] = stiffness3;
}

//-----------------------------------------------------------------------------

StiffnessQuad::StiffnessQuad(double stiffness1,
                             double stiffness2,
                             double stiffness3,
                             double stiffness4)
: _stiffness(std::vector<double>(MAX_NUM_STIFFNESS_THREADS, 0.0))
{
  _stiffness[0] = stiffness1;
  _stiffness[1] = stiffness2;
  _stiffness[2] = stiffness3;
  _stiffness[3] = stiffness4;
}

//-----------------------------------------------------------------------------

double StiffnessQuad::getStiffness(const int index) const {
  if (index < 0 || index >= MAX_NUM_STIFFNESS_THREADS)
    throw IndexAccessException();
  else return _stiffness[index];
}

//-----------------------------------------------------------------------------

void StiffnessQuad::setStiffness(const int index, const double value) {
  if (index < 0 || index >= MAX_NUM_STIFFNESS_THREADS)
    throw IndexAccessException();
  else _stiffness[index] = value;
}

//-----------------------------------------------------------------------------

void StiffnessQuad::setStiffness(const double stiffness1,
                                 const double stiffness2,
                                 const double stiffness3,
                                 const double stiffness4)
{
  _stiffness[0] = stiffness1;
  _stiffness[1] = stiffness2;
  _stiffness[2] = stiffness3;
  _stiffness[3] = stiffness4;
}

//-----------------------------------------------------------------------------
// PointInfo implementation
//-----------------------------------------------------------------------------

PointInfo::PointInfo()
: _load(LoadInfo()),
  _x(0.0),
  _y(0.0),
  _z(0.0),
  _fixed(false)
{}

//-----------------------------------------------------------------------------

PointInfo::PointInfo(double xArg, double yArg, double zArg, bool fixedArg)
: _load(LoadInfo()),
  _x(xArg),
  _y(yArg),
  _z(zArg),
  _fixed(fixedArg)
{}

//-----------------------------------------------------------------------------

PointInfo::PointInfo(double xArg,
                     double yArg,
                     double zArg,
                     const LoadInfo &loadArg,
                     bool fixedArg)
: _load(loadArg),
  _x(xArg),
  _y(yArg),
  _z(zArg),
  _fixed(fixedArg)
{}

//-----------------------------------------------------------------------------

void PointInfo::setPosition(const double xArg,
                            const double yArg,
                            const double zArg)
{
  _x = xArg;
  _y = yArg;
  _z = zArg;
}

//-----------------------------------------------------------------------------

void PointInfo::getPosition(double &xArg,
                            double &yArg,
                            double &zArg) const
{
  xArg = _x;
  yArg = _y;
  zArg = _z;
}

//-----------------------------------------------------------------------------

double PointInfo::x() const {
  return _x;
}

//-----------------------------------------------------------------------------

double PointInfo::y() const {
  return _y;
}

//-----------------------------------------------------------------------------

double PointInfo::z() const {
  return _z;
}

//-----------------------------------------------------------------------------

void PointInfo::setX(double value) {
  _x = value;
}

//-----------------------------------------------------------------------------

void PointInfo::setY(double value) {
  _y = value;
}

//-----------------------------------------------------------------------------

void PointInfo::setZ(double value) {
  _z = value;
}

//-----------------------------------------------------------------------------

void PointInfo::setFixed(bool fixArg) {
  _fixed = fixArg;
}

//-----------------------------------------------------------------------------

bool PointInfo::getFixed() const {
  return _fixed;
}

//-----------------------------------------------------------------------------

void PointInfo::setLoadInfo(const LoadInfo &loadArg) {
  _load = loadArg;
}

//-----------------------------------------------------------------------------

void PointInfo::setLoad(const double x, const double y, const double z) {
  _load.setLoad(x, y, z);
}

//-----------------------------------------------------------------------------

const LoadInfo &PointInfo::getLoadInfo() const {
  return _load;
}

//-----------------------------------------------------------------------------

LoadInfo &PointInfo::getLoadInfo() {
  return _load;
}

//-----------------------------------------------------------------------------
// EdgeInfo implementation
//-----------------------------------------------------------------------------

EdgeInfo::EdgeInfo()
: _ib(-1),
  _ie(-1),
  _stiffness(StiffnessQuad())
{}

//-----------------------------------------------------------------------------

EdgeInfo::EdgeInfo(const int ibArg,
                   const int ieArg,
                   const StiffnessQuad &stiffArg)
: _ib(ibArg),
  _ie(ieArg),
  _stiffness(stiffArg)
{}

//-----------------------------------------------------------------------------

const int EdgeInfo::ib() const {
  return _ib;
}

//-----------------------------------------------------------------------------

const int EdgeInfo::ie() const {
  return _ie;
}

//-----------------------------------------------------------------------------

void EdgeInfo::setIB(const int value) {
  _ib = value;
}

//-----------------------------------------------------------------------------

void EdgeInfo::setIE(const int value) {
  _ie = value;
}

//-----------------------------------------------------------------------------

void EdgeInfo::setStiffnessQuad(const StiffnessQuad &quad) {
  _stiffness = quad;
}

//-----------------------------------------------------------------------------

StiffnessQuad &EdgeInfo::getStiffnessQuad() {
  return _stiffness;
}

//-----------------------------------------------------------------------------

const StiffnessQuad &EdgeInfo::getStiffnessQuad() const {
  return _stiffness;
}

//-----------------------------------------------------------------------------

void EdgeInfo::setStiffness(const int index, const double value) {
  _stiffness.setStiffness(index, value);
}

//-----------------------------------------------------------------------------

void EdgeInfo::setStiffness(const double stiffness1,
                            const double stiffness2,
                            const double stiffness3,
                            const double stiffness4)
{
  _stiffness.setStiffness(
    stiffness1, stiffness2, stiffness3, stiffness4);
}

//-----------------------------------------------------------------------------

double EdgeInfo::getStiffness(const int index) const {
  return _stiffness.getStiffness(index);
}

//-----------------------------------------------------------------------------
// ForceDensity implementation
//-----------------------------------------------------------------------------

ForceDensity::ForceDensity(const int numThreads)
: _numThreads(std::max<const int>(std::min<const int>(
     numThreads, MAX_NUM_STIFFNESS_THREADS), 0)),
  _enabledInstance(true),
  _enabledThreads(std::bitset<MAX_NUM_STIFFNESS_THREADS>(0xF))
{
  for (int i = 0; i < _numThreads; ++i) {
    _loads.push_back(LoadInfo());
    _result.push_back(Eigen::MatrixXd());  
  }
}

//-----------------------------------------------------------------------------

int ForceDensity::getNumStiffnessThreads() const {
  return _numThreads;
}

//-----------------------------------------------------------------------------

bool ForceDensity::isEnabledInstance() const {
  return _enabledInstance;
}

//-----------------------------------------------------------------------------

void ForceDensity::enableInstance(const bool enable) {
  _enabledInstance = enable;
}

//-----------------------------------------------------------------------------

void ForceDensity::enableStiffnessThread(const int threadID,
                                         const bool enable)
{
  if (threadID >= 0 && threadID < MAX_NUM_THREADS)
    _enabledThreads.set(threadID, enable);
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

const ForceDensity::PointInfoVectorTy &ForceDensity::getPoints() const {
  return _points;
}

//-----------------------------------------------------------------------------

ForceDensity::PointInfoVectorTy &ForceDensity::getPoints() {
  return _points;
}

//-----------------------------------------------------------------------------

const ForceDensity::EdgeInfoVectorTy &ForceDensity::getEdges() const {
  return _edges;
}

//-----------------------------------------------------------------------------

ForceDensity::EdgeInfoVectorTy &ForceDensity::getEdges() {
  return _edges;
}

//-----------------------------------------------------------------------------

const ForceDensity::LoadVectorTy &ForceDensity::getUniformLoads() const {
  return _loads;
}

//-----------------------------------------------------------------------------

ForceDensity::LoadVectorTy &ForceDensity::getUniformLoads() {
  return _loads;
}

//-----------------------------------------------------------------------------

const ForceDensity::ResultVectorTy &ForceDensity::getResults() const {
  return _result;
}

//-----------------------------------------------------------------------------

ForceDensity::ResultVectorTy &ForceDensity::getResults() {
  return _result;
}

//-----------------------------------------------------------------------------

void ForceDensity::addPoint(const double x,
                            const double y,
                            const double z,
                            const bool fixed)
{
  _points.push_back(PointInfo(x, y, z, fixed));
}

//-----------------------------------------------------------------------------

void ForceDensity::addPoint(const double x,
                            const double y,
                            const double z,
                            const LoadInfo &load,
                            const bool fixed)
{
  _points.push_back(PointInfo(x, y, z, load, fixed));
}

//-----------------------------------------------------------------------------

PointInfo &ForceDensity::getPoint(const int index) {
  if (index >= 0 && index < _points.size())
    return _points[index];
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

const PointInfo &ForceDensity::getPoint(const int index) const {
  if (index >= 0 && index < _points.size())
    return _points[index];
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

void ForceDensity::setPoint(const int index, const PointInfo &p) {
  if (index >= 0 && index < _points.size())
    _points[index] = p;
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

void ForceDensity::addEdge(const int ib, 
                           const int ie,
                           const StiffnessQuad &stiffness)
{
  _edges.push_back(EdgeInfo(ib, ie, stiffness));
}

//-----------------------------------------------------------------------------

EdgeInfo &ForceDensity::getEdge(const int index) {
  if (index >= 0 && index < _edges.size())
    return _edges[index];
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

const EdgeInfo &ForceDensity::getEdge(const int index) const {
  if (index >= 0 && index < _edges.size())
    return _edges[index];
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

void ForceDensity::setEdge(const int index, const EdgeInfo &e) {
  if (index >= 0 && index < _edges.size())
    _edges[index] = e;
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

const LoadInfo &ForceDensity::getUniformLoad(const int index) const {
  if (index >= 0 && index < MAX_NUM_STIFFNESS_THREADS)
    return _loads[index];
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

LoadInfo &ForceDensity::getUniformLoad(const int index) {
  if (index >= 0 && index < MAX_NUM_STIFFNESS_THREADS)
    return _loads[index];
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

void ForceDensity::setUniformLoad(const int index, const LoadInfo &load) {
  if (index >= 0 && index < MAX_NUM_STIFFNESS_THREADS)
    _loads[index] = load;
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

void ForceDensity::setUniformLoad(const int index,
                                  const double loadX,
                                  const double loadY,
                                  const double loadZ)
{
  if (index >= 0 && index < MAX_NUM_STIFFNESS_THREADS)
    _loads[index].setLoad(loadX, loadY, loadZ);
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

const Eigen::MatrixXd &ForceDensity::getResult(const int index) const {
  if (index >= 0 && index < _result.size())
    return _result[index];
  else throw IndexAccessException();
}

//-----------------------------------------------------------------------------

void ForceDensity::clearInput() {
  _points.clear();
  _edges.clear();
}

//-----------------------------------------------------------------------------

void ForceDensity::setInput(const int threadID) {
  if (threadID < 0 || threadID >= _numThreads)
    throw IndexAccessException();

  for (int i = 0; i < _points.size(); ++i) {
    _points[i].setPosition(_result[threadID](i, 0),
                           _result[threadID](i, 1),
                           _result[threadID](i, 2));
  }
}

//-----------------------------------------------------------------------------

void ForceDensity::calculate() {
  std::list<std::thread> threads;

  for (int threadID = 0; threadID < _numThreads; ++threadID) {
    if (_enabledThreads.test(threadID) == true) {
      threads.push_back(
        std::thread(ForceDensityThread(),
                    threadID,
                    std::cref(_points),
                    std::cref(_edges),
                    std::cref(_loads[threadID]),
                    std::ref(_result[threadID])));
    }
  }

  for (auto &t : threads) t.join();
}

//-----------------------------------------------------------------------------
// ForceDensityThread implementation
//-----------------------------------------------------------------------------

void ForceDensityThread::operator()(
                       const int threadID,
                       const ForceDensity::PointInfoVectorTy &points,
                       const ForceDensity::EdgeInfoVectorTy &edges,
                       const LoadInfo &uniLoad,
                       Eigen::MatrixXd &x)
{
  const unsigned numEdges = edges.size();
  const unsigned numVertices = points.size();

  if (!numEdges || !numVertices) throw BadInputException();

  Eigen::SparseMatrix<double> C_free(numEdges, numVertices);
  Eigen::MatrixXd C_fixed =
    Eigen::MatrixXd::Zero(numEdges, numVertices);

  for (int E = 0; E < edges.size(); ++E) {
    if (edges[E].ib() == edges[E].ie()) throw BadEdgeException();
    if (edges[E].ib() >= numVertices) throw IndexAccessException();
    if (edges[E].ie() >= numVertices) throw IndexAccessException();

    if (!points[edges[E].ib()].getFixed())
      C_free.insert(E, edges[E].ib()) = -1;
    else C_fixed(E, edges[E].ib()) = -1;

    if (!points[edges[E].ie()].getFixed())
      C_free.insert(E, edges[E].ie()) = 1;
    else C_fixed(E, edges[E].ie()) = 1;
  }

  Eigen::SparseMatrix<double> Q(numEdges, numEdges);

  for (int E = 0; E < edges.size(); ++E) {
    Q.insert(E, E) = edges[E].getStiffness(threadID);
  }

  const Eigen::SparseMatrix<double> C_free_t = C_free.transpose();
  const Eigen::SparseMatrix<double> A = C_free_t * Q * C_free;

  Eigen::MatrixXd xyz_fixed = Eigen::MatrixXd::Zero(numVertices, 3);
  Eigen::MatrixXd p = Eigen::MatrixXd::Zero(numVertices, 3);

  std::list<unsigned> cornerIDs;

  for (unsigned V = 0; V < numVertices; ++V) {
    p(V, 0) = uniLoad.loadX() + points[V].getLoadInfo().loadX();
    p(V, 1) = uniLoad.loadY() + points[V].getLoadInfo().loadY();
    p(V, 2) = uniLoad.loadZ() + points[V].getLoadInfo().loadZ();

    if (points[V].getFixed()) {
      const PointInfo &corner = points[V];

      cornerIDs.push_back(V);

      xyz_fixed(V, 0) = corner.x();
      xyz_fixed(V, 1) = corner.y();
      xyz_fixed(V, 2) = corner.z();
    }
  }

  const Eigen::MatrixXd b = p - (C_free_t * Q * C_fixed * xyz_fixed);

  // now create a solver instance and run it. For the time being we stick to
  // a QR solver since no 'special' things are required to be met. A cholesky
  // solver might be a better choice in future, which however requires the 'A'
  // matrix to be sparse and positive definite, which in turn probably requi-
  // res a different build-up of A. This said, as alternative an QR solver
  // works just fine for now
  Eigen::SparseQR<
    Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > QR;

  QR.compute(A);

  if (QR.info() != Eigen::Success) throw SolverComputeException();

  x = QR.solve(b);

  if (QR.info() != Eigen::Success) throw SolverException();

  // for completeness, transfer the fixed points (as contained in the input)
  // array to the result matrix
  for (unsigned cornerID : cornerIDs) {
    x(cornerID, 0) = points[cornerID].x();
    x(cornerID, 1) = points[cornerID].y();
    x(cornerID, 2) = points[cornerID].z();
  }
}

//-----------------------------------------------------------------------------
// ForceDensityInstance implementation
//-----------------------------------------------------------------------------

void ForceDensityInstance::operator()(ForceDensity &fd) {
  fd.calculate();
}

//-----------------------------------------------------------------------------
// ForceDensityPool implementation
//-----------------------------------------------------------------------------

ForceDensityPool::ForceDensityPool()
{}

//-----------------------------------------------------------------------------

std::vector<ForceDensity*> &ForceDensityPool::getInstances() {
  return _instances;
}

//-----------------------------------------------------------------------------

const std::vector<ForceDensity*> &ForceDensityPool::getInstances() const {
  return _instances;
}

//-----------------------------------------------------------------------------

void ForceDensityPool::createInstance(const int numStiffnessThreads) {
  _instances.push_back(new ForceDensity(numStiffnessThreads));
}

//-----------------------------------------------------------------------------

void ForceDensityPool::addInstance(ForceDensity *instance) {
  if (instance == 0) throw BadObjectException();

  _instances.push_back(instance);
}

//-----------------------------------------------------------------------------

ForceDensity *ForceDensityPool::getInstance(const int instanceID) {
  if (instanceID < 0 || instanceID >= _instances.size())
    throw IndexAccessException();

  return _instances[instanceID];
}

//-----------------------------------------------------------------------------

void ForceDensityPool::releaseInstance(const int instanceID) {
  if (instanceID < 0 || instanceID >= _instances.size())
    throw IndexAccessException();

  ForceDensity *fd = _instances[instanceID];
  _instances[instanceID] = 0;

  delete fd;
}

//-----------------------------------------------------------------------------

void ForceDensityPool::releaseInstances() {
  for (auto *instance: _instances)
    delete instance;

  _instances.clear();
}

//-----------------------------------------------------------------------------

int ForceDensityPool::getNumInstances() {
  return _instances.size();
}

//-----------------------------------------------------------------------------

void ForceDensityPool::addPoint(const int instanceID,
                                const double x,
                                const double y,
                                const double z,
                                const bool fixed)
{
  if (instanceID < 0 || instanceID >= _instances.size())
    throw IndexAccessException();

  if (!_instances[instanceID])
    throw BadObjectException();

  _instances[instanceID]->addPoint(x, y, z, fixed);
}

//-----------------------------------------------------------------------------

void ForceDensityPool::addEdge(const int instanceID,
                               const int ib,
                               const int ie,
                               const double stiffness)
{
  if (instanceID < 0 || instanceID >= _instances.size())
    throw IndexAccessException();

  if (!_instances[instanceID])
    throw BadObjectException();

  _instances[instanceID]->addEdge(ib, ie, stiffness);
}

//-----------------------------------------------------------------------------

PointInfo &ForceDensityPool::getPoint(const int instanceID,
                                      const int pointIndex)
{
  if (instanceID < 0 || instanceID >= _instances.size())
    throw IndexAccessException();

  if (!_instances[instanceID])
    throw BadObjectException();

  return _instances[instanceID]->getPoint(pointIndex);
}

//-----------------------------------------------------------------------------

const PointInfo &ForceDensityPool::getPoint(const int instanceID,
                                            const int pointIndex) const
{
  if (instanceID < 0 || instanceID >= _instances.size())
    throw IndexAccessException();

  if (!_instances[instanceID])
    throw BadObjectException();

  return _instances[instanceID]->getPoint(pointIndex);
}

//-----------------------------------------------------------------------------

EdgeInfo &ForceDensityPool::getEdge(const int instanceID,
                                    const int edgeIndex)
{
  if (instanceID < 0 || instanceID >= _instances.size())
    throw IndexAccessException();

  if (!_instances[instanceID])
    throw BadObjectException();

  return _instances[instanceID]->getEdge(edgeIndex);
}

//-----------------------------------------------------------------------------

const EdgeInfo &ForceDensityPool::getEdge(const int instanceID,
                                          const int edgeIndex) const
{
  if (instanceID < 0 || instanceID >= _instances.size())
    throw IndexAccessException();

  if (!_instances[instanceID])
    throw BadObjectException();

  return _instances[instanceID]->getEdge(edgeIndex);
}

//-----------------------------------------------------------------------------

void ForceDensityPool::runInstances() {
  std::list<std::thread> threads;

  for (ForceDensity *fd : _instances) {
    threads.push_back(
      std::thread(ForceDensityInstance(),
                  std::ref(*fd)));
  }

  for (auto &t : threads) t.join();
}

