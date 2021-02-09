#include "algorithm.h"
#include "forcedensity.h"

#include <iostream>
#include <thread>

#define DEBUG_PRINTS

/******************************************************************************/

int main(int argc, char **argv) {
  if (createForceDensityPool() != NO_ERROR) return 1;

  int mesh1, mesh2;

  // mesh 1 data, 1 stiffness thread
  {
    if (addForceDensityInstance(1, mesh1) != NO_ERROR) return 1;

    addInstancePointInfo(mesh1, -1.0, -1.0, 0.0, true);  // a: 0
    addInstancePointInfo(mesh1, -1.0, 0.0, 0.0, false);  // b: 1
    addInstancePointInfo(mesh1, -1.0, 1.0, 0.0, true);   // c: 2
    addInstancePointInfo(mesh1, 0.0, 1.0, 0.0, false);   // d: 3
    addInstancePointInfo(mesh1, 1.0, 1.0, 0.0, true);    // e: 4
    addInstancePointInfo(mesh1, 1.0, 0.0, 0.0, false);   // f: 5
    addInstancePointInfo(mesh1, 1.0, -1.0, 0.0, true);   // g: 6
    addInstancePointInfo(mesh1, 0.0, -1.0, 0.0, false);  // h: 7
    addInstancePointInfo(mesh1, 0.0, 0.0, 0.0, false);   // i: 8

    addInstanceEdgeInfo(mesh1, 0, 1, 1.0, 0, 0, 0);   // a->b
    addInstanceEdgeInfo(mesh1, 1, 2, 1.0, 0, 0, 0);   // b->c
    addInstanceEdgeInfo(mesh1, 2, 3, 1.0, 0, 0, 0);   // c->d
    addInstanceEdgeInfo(mesh1, 3, 4, 1.0, 0, 0, 0);   // d->e
    addInstanceEdgeInfo(mesh1, 4, 5, 1.0, 0, 0, 0);   // e->f
    addInstanceEdgeInfo(mesh1, 5, 6, 1.0, 0, 0, 0);   // f->g
    addInstanceEdgeInfo(mesh1, 6, 7, 10.0, 0, 0, 0);  // g->h
    addInstanceEdgeInfo(mesh1, 7, 0, 10.0, 0, 0, 0);  // h->a
    addInstanceEdgeInfo(mesh1, 1, 8, 10.0, 0, 0, 0);  // b->i
    addInstanceEdgeInfo(mesh1, 3, 8, 10.0, 0, 0, 0);  // d->i
    addInstanceEdgeInfo(mesh1, 5, 8, 10.0, 0, 0, 0);  // f->i
    addInstanceEdgeInfo(mesh1, 7, 8, 10.0, 0, 0, 0);  // h->i
  }

  // mesh 2 data, 4 stiffness threads
  {
    // add an instance for mesh 2 (ID 1) with 4 stiffness values
    if (addForceDensityInstance(4, mesh2) != NO_ERROR) return 1;

    addInstancePointInfo(mesh2, -1.0, -1.0, 0.0, true);  // a: 0
    addInstancePointInfo(mesh2, -1.0, 0.0, 0.0, false);  // b: 1
    addInstancePointInfo(mesh2, -1.0, 1.0, 0.0, true);   // c: 2
    addInstancePointInfo(mesh2, 0.0, 1.0, 0.0, false);   // d: 3
    addInstancePointInfo(mesh2, 1.0, 1.0, 0.0, true);    // e: 4
    addInstancePointInfo(mesh2, 1.0, 0.0, 0.0, false);   // f: 5
    addInstancePointInfo(mesh2, 1.0, -1.0, 0.0, true);   // g: 6
    addInstancePointInfo(mesh2, 0.0, -1.0, 0.0, false);  // h: 7
    addInstancePointInfo(mesh2, 0.0, 0.0, 0.0, false);   // i: 8
  
    addInstanceEdgeInfo(mesh2, 0, 1, 5.0, 10.0, 12.0, 2.0);  // a->b
    addInstanceEdgeInfo(mesh2, 1, 2, 5.0, 10.0, 12.0, 2.0);  // b->c
    addInstanceEdgeInfo(mesh2, 2, 3, 5.0, 10.0, 12.0, 2.0);  // c->d
    addInstanceEdgeInfo(mesh2, 3, 4, 5.0, 10.0, 12.0, 2.0);  // d->e
    addInstanceEdgeInfo(mesh2, 4, 5, 5.0, 10.0, 12.0, 2.0);  // e->f
    addInstanceEdgeInfo(mesh2, 5, 6, 5.0, 10.0, 12.0, 2.0);  // f->g
    addInstanceEdgeInfo(mesh2, 6, 7, 5.0, 10.0, 12.0, 2.0);  // g->h
    addInstanceEdgeInfo(mesh2, 7, 0, 5.0, 10.0, 12.0, 2.0);  // h->a
    addInstanceEdgeInfo(mesh2, 1, 8, 5.0, 10.0, 12.0, 2.0);  // b->i
    addInstanceEdgeInfo(mesh2, 3, 8, 5.0, 10.0, 12.0, 2.0);  // d->i
    addInstanceEdgeInfo(mesh2, 5, 8, 5.0, 10.0, 12.0, 2.0);  // f->i
    addInstanceEdgeInfo(mesh2, 7, 8, 5.0, 10.0, 12.0, 2.0);  // h->i
  }

  // enable/disable particular stiffness threads, all 4 threads are
  // initially enabled by default
  enableInstanceThread(mesh2, 1, false);
  enableInstanceThread(mesh2, 3, false);

  // supply uniform surface loads per instance/thread, thus, up to 4
  // load conditions can be computed per instance at once
  setInstanceUniformLoad(mesh2, 0, 10.0, 10.0, 1.0);
  setInstanceUniformLoad(mesh2, 2, -1.0, 5.0, 10.0);

  // execute all instances (meshes) in parallel. Mesh 2 additionally
  // triggers 4 threads, one for each edge/stiffness variant
  if (runForceDensityInstances() != NO_ERROR) return 1;

  // result of mesh 1 and stiffness thread 0
  {
    for (int i = 0; i < getNumInstancePoints(mesh1); ++i) {
      double x, y, z;

      // only one stiffness-thread has been specified, thus accessing
      // any other threads will result in an error
      getInstanceResultPoint(mesh1, 0, i, x, y, z);
    }
  }

  // result of mesh 2 and stiffness threads 0 and 2. Threads 1 and 3
  // have been disabled, see above
  {
    for (int i = 0; i < getNumInstancePoints(mesh2); ++i) {
      double x_t0, y_t0, z_t0;
      double x_t2, y_t2, z_t2;

      getInstanceResultPoint(mesh2, 0, i, x_t0, y_t0, z_t0);
      getInstanceResultPoint(mesh2, 2, i, x_t2, y_t2, z_t2);
    }
  }

#ifdef DEBUG_PRINTS
  // visualize output for thread 1 of mesh 1
  std::cout << "  result of thread 0 for mesh "
            << mesh1 << ":" << std::endl;

  for (int i = 0; i < getNumInstancePoints(mesh1); ++i) {
    double x, y, z;
    getInstanceResultPoint(mesh1, 0, i, x, y, z);
    
    std::cout << "    point " << i << ": ("
              << x << ", " << y << ", " << z << ")"
              << std::endl;
  }

  // visualize output for thread 0 of mesh 2
  std::cout << std::endl << "  result of thread 0 for mesh "
            << mesh2 << ":" << std::endl;

  for (int i = 0; i < getNumInstancePoints(mesh2); ++i) {
    double x, y, z;
    getInstanceResultPoint(mesh2, 0, i, x, y, z);

    std::cout << "    point " << i << ": ("
              << x << ", " << y << ", " << z << ")"
              << std::endl;
  }

  // visualize output for thread 2 of mesh 2
  std::cout << std::endl << "  result of thread 2 for mesh "
            << mesh2 << ":" << std::endl;

  for (int i = 0; i < getNumInstancePoints(mesh2); ++i) {
    double x, y, z;
    getInstanceResultPoint(mesh2, 2, i, x, y, z);

    std::cout << "    point " << i << ": ("
              << x << ", " << y << ", " << z << ")"
              << std::endl;
  }
#endif

  // re-enable all threads for mesh 2 
  enableInstanceThread(mesh2, 1, true);
  enableInstanceThread(mesh2, 3, true);

  // skip mesh 1 computation, this has no runtime overhead,
  // but occupies memory (until the instance is dropped)
  enableForceDensityInstance(mesh1, false);

  // run again, mesh 2 triggers all 4 stiffness threads now
  if (runForceDensityInstances() != NO_ERROR) return 1;

  // done, goodbye
  releaseForceDensityPool();
  return 0;
}
