#ifndef ALGORITHM_H
#define ALGORITHM_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <vector>
#include <bitset>

const int MAX_NUM_STIFFNESS_THREADS = 4;

//-----------------------------------------------------------------------------
// @enum ErrorCode
// @brief Capture error cases expected to happen most during execution

enum ErrorCode {
  NO_ERROR = 0,

  // general error not yet associated with any special case, this is used
  // only once and actually should never happen
  UNKNOWN_ERROR,

  // for the time being (simplicity) allow only one thread pool at a time
  SINGLE_POOL_ERROR,

  // happens when trying to access an instance that does not yet exists
  BAD_INSTANCE_ERROR,

  // point/edge index out of bounds
  BAD_INDEX_ERROR,

  // specificially provided to catch edges with start and end indices being
  // equal
  MALFORMED_EDGE_ERROR,

  // no data (empty point/edge containers) specified for the solver
  EMPTY_ALGORITHM_INPUT_ERROR,

  // point/edge index out of bounds within the algorithm itself. Should
  // actually never happen, but i leave it in place for debugging yet
  BAD_ALGORITHM_INDEX_ERROR,

  // aggregate error for 'Eigen::compute' prior to equation solving
  BAD_EIGEN_COMPUTE_ERROR,

  // happens when the computation of the coefficient matrix A succeeds,
  // but the final solver process fails due to unexpected errors
  BAD_EIGEN_SOLVE_ERROR
};

//-----------------------------------------------------------------------------
// @class LoadInfo
// @brief Simple triple of doubles to capture load values

class LoadInfo {
  protected:
    double _loadX;
    double _loadY;
    double _loadZ;

  public:
    LoadInfo();
    LoadInfo(double x, double y, double z);

    double loadX() const;
    double loadY() const;
    double loadZ() const;

    void setLoadX(double);
    void setLoadY(double);
    void setLoadZ(double);

    void getLoad(double &x, double &y, double &z);
    void setLoad(double x, double y, double z);
};

//-----------------------------------------------------------------------------
// @class StiffnessQuad 
// @brief Simple quad of doubles to capture max. count of stiffness values

class StiffnessQuad {
  protected:
    std::vector<double> _stiffness;

  public:
    StiffnessQuad();
    StiffnessQuad(double stiffness);

    StiffnessQuad(double stiffness1,
                  double stiffness2);

    StiffnessQuad(double stiffness1,
                  double stiffness2,
                  double stiffness3);

    StiffnessQuad(double stiffness1,
                  double stiffness2,
                  double stiffness3,
                  double stiffness4);

    double getStiffness(const int index) const;
    void setStiffness(const int index, const double value);
    void setStiffness(const double stiffness1,
                      const double stiffness2,
                      const double stiffness3,
                      const double stiffness4);
};

//-----------------------------------------------------------------------------
// @class PointInfo
// @brief Capture most necessary information for a given point

class PointInfo {
     LoadInfo _load;

     double _x;
     double _y;
     double _z;
     bool _fixed;

  public:
    PointInfo();

    PointInfo(double x,
              double y,
              double z,
              bool fixed=false);

    PointInfo(double x,
              double y,
              double z,
              const LoadInfo &load,
              bool fixed=false);

    void getPosition(double &x, double &y, double &z) const;
    void setPosition(const double x, const double y, const double z);

    void setFixed(bool);
    bool getFixed() const;

    double x() const;
    double y() const;
    double z() const;

    void setX(double);
    void setY(double);
    void setZ(double);

    void setLoadInfo(const LoadInfo &load);

    void setLoad(const double loadX,
                 const double loadY,
                 const double loadZ);

    const LoadInfo &getLoadInfo() const;
    LoadInfo &getLoadInfo();
};

//-----------------------------------------------------------------------------
// @class EdgeInfo
// @brief Represent a single edge as a pair of vertex indices

class EdgeInfo {
    // begin and end index of the points this edge is incident to. Despite the
    // name which is an artifact we dont yet consider edges being directed
    int _ib;
    int _ie;

    // for additional shapes associate a stiffness value with each edge. This
    // value is usually relevant for outline edges and is 1 for points within
    // the surface
    StiffnessQuad _stiffness;

  public:
    EdgeInfo();
    EdgeInfo(const int ibArg, const int ieArg, const StiffnessQuad &stiff);

    // read only successors for the point indices for consistency. Thus, these
    // are meant to be set by the ctor and remain unchanged through the entire
    // lifespan of the info
    const int ib() const;
    const int ie() const;

    void setIB(const int);
    void setIE(const int);

    void setStiffnessQuad(const StiffnessQuad &quad);
    StiffnessQuad &getStiffnessQuad();
    const StiffnessQuad &getStiffnessQuad() const;

    double getStiffness(const int index) const;
    void setStiffness(const int index, const double value);
    void setStiffness(const double stiffness1,
                      const double stiffness2,
                      const double stiffness3,
                      const double stiffness4);
};

//-----------------------------------------------------------------------------
// @class ForceDensity
// @brief A simple variant of the force density algorithm for multithreaded
//   fixity calculation

class ForceDensity {
  public:
    static const int MAX_NUM_THREADS = 4;

    typedef std::vector<PointInfo> PointInfoVectorTy;
    typedef std::vector<EdgeInfo> EdgeInfoVectorTy;
    typedef std::vector<Eigen::MatrixXd> ResultVectorTy;
    typedef std::vector<LoadInfo> LoadVectorTy;
    typedef std::bitset<MAX_NUM_THREADS> EnableBitsetTy;

  protected:
    const int _numThreads;
    bool _enabledInstance;

    PointInfoVectorTy _points;
    EdgeInfoVectorTy _edges;
    LoadVectorTy _loads;
    ResultVectorTy _result;
    EnableBitsetTy _enabledThreads;

  public:
    ForceDensity(const int numThreads);

    int getNumStiffnessThreads() const;

    bool isEnabledInstance() const;
    void enableInstance(const bool enable);

    const PointInfoVectorTy &getPoints() const;
    PointInfoVectorTy &getPoints();

    const EdgeInfoVectorTy &getEdges() const;
    EdgeInfoVectorTy &getEdges();

    const LoadVectorTy &getUniformLoads() const;
    LoadVectorTy &getUniformLoads();

    const ResultVectorTy &getResults() const;
    ResultVectorTy &getResults();

    PointInfo &getPoint(const int index);
    const PointInfo &getPoint(const int index) const;
    void setPoint(const int index, const PointInfo &p);

    EdgeInfo &getEdge(const int index);
    const EdgeInfo &getEdge(const int index) const;
    void setEdge(const int index, const EdgeInfo &e);

    const LoadInfo &getUniformLoad(const int index) const;
    LoadInfo &getUniformLoad(const int index);
    void setUniformLoad(const int index, const LoadInfo &load);
    void setUniformLoad(const int index,
                        const double loadX,
                        const double loadY,
                        const double loadZ);

    const Eigen::MatrixXd &getResult(const int index) const;

    void addPoint(const double x,
                  const double y,
                  const double z,
                  const bool fixed = false);

    void addPoint(const double x,
                  const double y,
                  const double z,
                  const LoadInfo &load,
                  const bool fixed = false);

    void addEdge(const int ib,
                 const int ie,
                 const StiffnessQuad &stiffness);


    void enableStiffnessThread(const int threadID,
                               const bool enable);

    void clearInput();
    void setInput(const int threadID);
    void calculate();
};

//-----------------------------------------------------------------------------
// @class ForceDensityThread
// @brief A simple variant of the force density algorithm for form-finding

class ForceDensityThread {
  public:
    void operator()(const int threadID,
                    const ForceDensity::PointInfoVectorTy &points,
                    const ForceDensity::EdgeInfoVectorTy &edges,
                    const LoadInfo &uniLoad,
                    Eigen::MatrixXd &x);
};

//-----------------------------------------------------------------------------
// @class ForceDensityInstance
// @brief A simple functor for the top level threads

class ForceDensityInstance {
  public:
    void operator()(ForceDensity &fd);
};

//-----------------------------------------------------------------------------
// @class ForceDensityPool
// @brief A pool of high level force-density threads

class ForceDensityPool {
  protected:
    std::vector<ForceDensity*> _instances;

  public:
    ForceDensityPool();

    std::vector<ForceDensity*> &getInstances();
    const std::vector<ForceDensity*> &getInstances() const;

    void createInstance(const int numStiffnessThreads);
    void addInstance(ForceDensity *instance);

    void releaseInstance(const int instanceID);
    void releaseInstances();

    int getNumInstances();

    ForceDensity *getInstance(const int instanceID);

    // convenience interface for a particular instance
    void addPoint(const int instanceID,
                  const double x,
                  const double y,
                  const double z,
                  const bool fixed = false);

    void addEdge(const int instanceID,
                 const int ib,
                 const int ie,
                 const double stiffness = 10.0);

    PointInfo &getPoint(const int instanceID,
                        const int pointIndex);

    const PointInfo &getPoint(const int instanceID,
                              const int pointIndex) const;

    EdgeInfo &getEdge(const int instanceID,
                      const int edgeIndex);

    const EdgeInfo &getEdge(const int instanceID,
                            const int edgeIndex) const;

    void runInstances();
};

#endif // ALGORITHM_H
