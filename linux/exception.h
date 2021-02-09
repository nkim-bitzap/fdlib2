#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <exception>

//-----------------------------------------------------------------------------
// @enum ForceDensityException
// @brief Base class for exceptions within formfinding machinery 

class ForceDensityException : public std::exception {
  virtual const char *what() const throw() {
    return "Exception: general force-density error";
  }
};

//-----------------------------------------------------------------------------
// @enum BadObjectException 
// @brief Exception to catch bad (null) objects

class BadObjectException : public ForceDensityException {
  virtual const char *what() const throw() {
    return "Exception: invalid (null) object";
  }
};

//-----------------------------------------------------------------------------
// @enum IndexAccessException 
// @brief Exception to catch bad (array/container) accesses

class IndexAccessException : public ForceDensityException {
  virtual const char *what() const throw() {
    return "Exception: index access out of bounds";
  }
};

//-----------------------------------------------------------------------------
// @enum BadInputException
// @brief Exception to catch invalid (empty) point/edge data

class BadInputException : public ForceDensityException {
  virtual const char *what() const throw() {
    return "Exception: invalid algorithm input (points/edges)";
  }
};

//-----------------------------------------------------------------------------
// @enum BadEdgeException
// @brief Exception to catch circular edges (ie == ib)

class BadEdgeException : public ForceDensityException {
  virtual const char *what() const throw() {
    return "Exception: bad (circular) edge (ie == ib)";
  }
};

//-----------------------------------------------------------------------------
// @enum SolverComputeException
// @brief Exception to catch solver 'compute' state

class SolverComputeException : public ForceDensityException {
  virtual const char *what() const throw() {
    return "Exception: internal solver error (compute)";
  }
};

//-----------------------------------------------------------------------------
// @enum SolverException
// @brief Exception to catch solver 'solve' state

class SolverException : public ForceDensityException {
  virtual const char *what() const throw() {
    return "Exception: internal solver error (solve)";
  }
};

#endif /* EXCEPTION_H */

