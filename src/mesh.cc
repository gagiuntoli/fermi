#include "mesh.h"

#include "element.h"
#include "shape.h"

Mesh Mesh::create1Dlinear(size_t nx, double lx) {
  std::vector<Node> nodes;
  double h = lx / (nx - 1);
  for (size_t i = 0; i < nx; i++) {
    nodes.push_back({i * h, 0, 0});
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < nx - 1; i++) {
    std::vector<size_t> nodeIndexes = {i, i + 1};
    std::vector<Node> elemNodes;
    for (const size_t& index : nodeIndexes) {
      elemNodes.push_back(nodes[index]);
    }
    elements.push_back(std::make_shared<ElementDiffusion<1>>(std::make_shared<Segment2>(), elemNodes, nodeIndexes, 0.03,
                                                             0.04, 1.0, 1.2));
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements)};
}

Mesh Mesh::create2DlinearTria3(size_t nx, size_t ny, double lx, double ly) {
  std::vector<Node> nodes;
  double hx = lx / (nx - 1);
  double hy = ly / (ny - 1);
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      nodes.push_back({i * hx, j * hy, 0});
    }
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < (nx - 1); i++) {
    for (size_t j = 0; j < (ny - 1); j++) {
      std::vector<size_t> nodeIndexes = {i * ny + j, (i + 1) * ny + j, (i + 1) * ny + j + 1};
      std::vector<Node> elemNodes;
      for (const size_t& index : nodeIndexes) {
        elemNodes.push_back(nodes[index]);
      }
      elements.push_back(std::make_shared<ElementDiffusion<2>>(std::make_shared<Triangle3>(), elemNodes, nodeIndexes,
                                                               0.03, 0.04, 1.0, 1.2));

      nodeIndexes = {i * ny + j, i * ny + j + 1, (i + 1) * ny + j + 1};
      elemNodes.clear();
      for (const size_t& index : nodeIndexes) {
        elemNodes.push_back(nodes[index]);
      }
      elements.push_back(std::make_shared<ElementDiffusion<2>>(std::make_shared<Triangle3>(), elemNodes, nodeIndexes,
                                                               0.03, 0.04, 1.0, 1.2));
    }
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements)};
}

Mesh Mesh::create2DlinearQuad4(size_t nx, size_t ny, double lx, double ly) {
  std::vector<Node> nodes;
  double hx = lx / (nx - 1);
  double hy = ly / (ny - 1);
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      nodes.push_back({i * hx, j * hy, 0});
    }
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < (nx - 1); i++) {
    for (size_t j = 0; j < (ny - 1); j++) {
      std::vector<size_t> nodeIndexes = {i * ny + j, i * ny + j + 1, (i + 1) * ny + j, (i + 1) * ny + j + 1};
      std::vector<Node> elemNodes;
      for (const size_t& index : nodeIndexes) {
        elemNodes.push_back(nodes[index]);
      }
      elements.push_back(std::make_shared<ElementDiffusion<2>>(std::make_shared<Quadrilateral4>(), elemNodes,
                                                               nodeIndexes, 0.03, 0.04, 1.0, 1.2));
    }
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements)};
}

Mesh Mesh::create3DlinearHexa8(size_t nx, size_t ny, size_t nz, double lx, double ly, double lz) {
  std::vector<Node> nodes;
  double hx = lx / (nx - 1);
  double hy = ly / (ny - 1);
  double hz = lz / (nz - 1);
  for (size_t i = 0; i < nx; i++) {
    for (size_t j = 0; j < ny; j++) {
      for (size_t k = 0; k < nz; k++) {
        nodes.push_back({i * hx, j * hy, k * hz});
      }
    }
  }

  std::vector<std::shared_ptr<ElementBase>> elements;
  for (size_t i = 0; i < (nx - 1); i++) {
    for (size_t j = 0; j < (ny - 1); j++) {
      for (size_t k = 0; k < (nz - 1); k++) {
        std::vector<size_t> nodeIndexes = {i * ny + j + k * (nx * ny),
                                           i * ny + j + 1 + k * (nx * ny),
                                           (i + 1) * ny + j + k * (nx * ny),
                                           (i + 1) * ny + j + 1 + k * (nx * ny),  //
                                           i * ny + j + (k + 1) * (nx * ny),
                                           i * ny + j + 1 + (k + 1) * (nx * ny),
                                           (i + 1) * ny + j + (k + 1) * (nx * ny),
                                           (i + 1) * ny + j + 1 + (k + 1) * (nx * ny)};
        std::vector<Node> elemNodes;
        for (const size_t& index : nodeIndexes) {
          elemNodes.push_back(nodes[index]);
        }
        elements.push_back(std::make_shared<ElementDiffusion<3>>(std::make_shared<Hexagon8>(), elemNodes, nodeIndexes,
                                                                 0.03, 0.04, 1.0, 1.2));
      }
    }
  }

  return Mesh{.nodes = std::move(nodes), .elements = std::move(elements)};
}
