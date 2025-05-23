/*
 *  This source code is part of Fermi: a finite element code
 *  to solve the neutron diffusion problem for nuclear reactor
 *  designs.
 *
 *  Copyright (C) - 2019 - Guido Giuntoli <gagiuntoli@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published
 *  by the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef SHAPE_H
#define SHAPE_H

#include <array>
#include <iomanip>

#include "node.h"

inline static const size_t MAX_NGP = 8;
inline static const size_t MAX_SHAPES = 8;

typedef std::array<std::array<double, MAX_NGP>, MAX_SHAPES> Shapes;
typedef std::array<std::array<std::array<double, MAX_NGP>, 3>, MAX_SHAPES> Derivatives;

class Shape {
 public:
  constexpr Shape(size_t ngp_, const std::array<Node, MAX_NGP>& gps_, const std::array<double, MAX_NGP>& wgs_,
                  size_t nsh_, Shapes sh_, Derivatives dsh_)
      : ngp(ngp_), gps(gps_), wgs(wgs_), nsh(nsh_), sh(sh_), dsh(dsh_) {}
  const size_t ngp;
  const std::array<Node, MAX_NGP> gps;
  const std::array<double, MAX_NGP> wgs;
  const size_t nsh;
  const Shapes sh;
  const Derivatives dsh;

  std::string toString() const {
    std::ostringstream oss;
    oss << "gps:" << std::endl;
    for (size_t i = 0; i < ngp; i++) {
      oss << gps[i].toString() << std::endl;
    }

    oss << "wgs:" << std::endl;
    for (size_t i = 0; i < ngp; i++) {
      oss << wgs[i] << std::endl;
    }

    oss << "shapes:" << std::endl;
    oss << std::setw(5) << std::setprecision(5) << std::fixed;
    for (size_t i = 0; i < nsh; i++) {
      for (size_t gp = 0; gp < ngp; gp++) {
        oss << sh[i][gp] << (gp != ngp - 1 ? "," : "");
      }
      oss << std::endl;
    }

    oss << "dshapes" << std::endl;
    for (size_t i = 0; i < nsh; i++) {
      for (size_t d = 0; d < 3; d++) {
        oss << "sh: " << i << " d: " << d << " ";
        for (size_t gp = 0; gp < ngp; gp++) {
          oss << dsh[i][d][gp] << (gp != ngp - 1 ? "," : "");
        }
        oss << (d != 2 ? "\n" : "");
      }
      oss << (i != nsh - 1 ? "\n" : "");
    }
    return oss.str();
  }
};

class Segment2 : public Shape {
 public:
  constexpr Segment2() : Shape(2, makeGPS(), makeWGS(), 2, makeShapes(), makeDerivatives()) {}

  std::string toString() const {
    std::ostringstream oss;
    oss << "Segment Linear" << std::endl;
    oss << Shape::toString();
    return oss.str();
  }

 private:
  constexpr std::array<Node, MAX_NGP> makeGPS() {
    std::array<Node, MAX_NGP> arr{};
    arr[0] = {-0.577350269189626, 0, 0};
    arr[1] = {+0.577350269189626, 0, 0};
    return arr;
  }

  constexpr std::array<double, MAX_NGP> makeWGS() {
    std::array<double, MAX_NGP> arr{};
    arr[0] = 1.0;
    arr[1] = 1.0;
    return arr;
  }

  constexpr Shapes makeShapes() {
    Shapes sh{};
    for (size_t gp = 0; gp < 2; gp++) {
      sh[0][gp] = +0.5 * (1.0 - gps[gp].x);
      sh[1][gp] = +0.5 * (1.0 + gps[gp].x);
    }
    return sh;
  }

  constexpr Derivatives makeDerivatives() {
    Derivatives ds{};
    for (size_t gp = 0; gp < 2; gp++) {
      ds[0][0][gp] = -0.5;
      ds[1][0][gp] = +0.5;
    }
    return ds;
  }
};

class Triangle3 : public Shape {
 public:
  constexpr Triangle3() : Shape(3, makeGPS(), makeWGS(), 3, makeShapes(), makeDerivatives()) {}

  std::string toString() const {
    std::ostringstream oss;
    oss << "Triangle Linear" << std::endl;
    oss << Shape::toString();
    return oss.str();
  }

 private:
  constexpr std::array<Node, MAX_NGP> makeGPS() {
    std::array<Node, MAX_NGP> arr{};
    arr[0] = {+0.166666667, +0.166666667, 0};
    arr[1] = {+0.666666667, +0.166666667, 0};
    arr[2] = {+0.166666667, +0.666666667, 0};
    return arr;
  }

  constexpr std::array<double, MAX_NGP> makeWGS() {
    std::array<double, MAX_NGP> arr{};
    arr[0] = 0.166666667;
    arr[1] = 0.166666667;
    arr[2] = 0.166666667;
    return arr;
  }

  constexpr Shapes makeShapes() {
    Shapes sh{};
    for (size_t gp = 0; gp < 3; gp++) {
      sh[0][gp] = 1.0 - gps[gp].x - gps[gp].y;
      sh[1][gp] = gps[gp].x;
      sh[2][gp] = gps[gp].y;
    }
    return sh;
  }

  constexpr Derivatives makeDerivatives() {
    Derivatives ds{};
    for (size_t gp = 0; gp < 3; gp++) {
      ds[0][0][gp] = -1.0;
      ds[0][1][gp] = -1.0;

      ds[1][0][gp] = +1.0;
      ds[1][1][gp] = +0.0;

      ds[2][0][gp] = +0.0;
      ds[2][1][gp] = +1.0;
    }
    return ds;
  }
};

class Quadrilateral4 : public Shape {
 public:
  constexpr Quadrilateral4() : Shape(4, makeGPS(), makeWGS(), 4, makeShapes(), makeDerivatives()) {}

  std::string toString() const {
    std::ostringstream oss;
    oss << "Quad Linear" << std::endl;
    oss << Shape::toString();
    return oss.str();
  }

 private:
  constexpr std::array<Node, MAX_NGP> makeGPS() {
    std::array<Node, MAX_NGP> arr{};
    arr[0] = {-0.577350269189626, -0.577350269189626, 0};
    arr[1] = {-0.577350269189626, +0.577350269189626, 0};
    arr[2] = {+0.577350269189626, -0.577350269189626, 0};
    arr[3] = {+0.577350269189626, +0.577350269189626, 0};
    return arr;
  }

  constexpr std::array<double, MAX_NGP> makeWGS() {
    std::array<double, MAX_NGP> arr{};
    for (size_t gp = 0; gp < 4; gp++) {
      arr[gp] = 1.0;
    }
    return arr;
  }

  constexpr Shapes makeShapes() {
    Shapes sh{};
    for (size_t gp = 0; gp < 4; gp++) {
      sh[0][gp] = (1.0 - gps[gp].x) * (1.0 - gps[gp].y) * 0.25;
      sh[1][gp] = (1.0 - gps[gp].x) * (1.0 + gps[gp].y) * 0.25;
      sh[2][gp] = (1.0 + gps[gp].x) * (1.0 - gps[gp].y) * 0.25;
      sh[3][gp] = (1.0 + gps[gp].x) * (1.0 + gps[gp].y) * 0.25;
    }
    return sh;
  }

  constexpr Derivatives makeDerivatives() {
    Derivatives ds{};
    for (size_t gp = 0; gp < 4; gp++) {
      ds[0][0][gp] = -1.0 * (1.0 - gps[gp].y) * 0.25;
      ds[0][1][gp] = -1.0 * (1.0 - gps[gp].x) * 0.25;

      ds[1][0][gp] = -1.0 * (1.0 + gps[gp].y) * 0.25;
      ds[1][1][gp] = +1.0 * (1.0 - gps[gp].x) * 0.25;

      ds[2][0][gp] = +1.0 * (1.0 - gps[gp].y) * 0.25;
      ds[2][1][gp] = -1.0 * (1.0 + gps[gp].x) * 0.25;

      ds[3][0][gp] = +1.0 * (1.0 + gps[gp].y) * 0.25;
      ds[3][1][gp] = +1.0 * (1.0 + gps[gp].x) * 0.25;
    }
    return ds;
  }
};

class Hexagon8 : public Shape {
 public:
  constexpr Hexagon8() : Shape(8, makeGPS(), makeWGS(), 8, makeShapes(), makeDerivatives()) {}

  std::string toString() const {
    std::ostringstream oss;
    oss << "Hexagon Linear" << std::endl;
    oss << Shape::toString();
    return oss.str();
  }

 private:
  constexpr std::array<Node, MAX_NGP> makeGPS() {
    std::array<Node, MAX_NGP> arr{};
    arr[0] = {-0.577350269189626, -0.577350269189626, -0.577350269189626};
    arr[1] = {-0.577350269189626, +0.577350269189626, -0.577350269189626};
    arr[2] = {+0.577350269189626, -0.577350269189626, -0.577350269189626};
    arr[3] = {+0.577350269189626, +0.577350269189626, -0.577350269189626};
    arr[4] = {-0.577350269189626, -0.577350269189626, +0.577350269189626};
    arr[5] = {-0.577350269189626, +0.577350269189626, +0.577350269189626};
    arr[6] = {+0.577350269189626, -0.577350269189626, +0.577350269189626};
    arr[7] = {+0.577350269189626, +0.577350269189626, +0.577350269189626};
    return arr;
  }

  constexpr std::array<double, MAX_NGP> makeWGS() {
    std::array<double, MAX_NGP> arr{};
    for (size_t gp = 0; gp < 8; gp++) {
      arr[gp] = 1.0;
    }
    return arr;
  }

  constexpr Shapes makeShapes() {
    Shapes sh{};
    for (size_t gp = 0; gp < 8; gp++) {
      sh[0][gp] = (1.0 - gps[gp].x) * (1.0 - gps[gp].y) * (1.0 - gps[gp].z) * 0.125;
      sh[1][gp] = (1.0 - gps[gp].x) * (1.0 + gps[gp].y) * (1.0 - gps[gp].z) * 0.125;
      sh[2][gp] = (1.0 + gps[gp].x) * (1.0 - gps[gp].y) * (1.0 - gps[gp].z) * 0.125;
      sh[3][gp] = (1.0 + gps[gp].x) * (1.0 + gps[gp].y) * (1.0 - gps[gp].z) * 0.125;
      sh[4][gp] = (1.0 - gps[gp].x) * (1.0 - gps[gp].y) * (1.0 + gps[gp].z) * 0.125;
      sh[5][gp] = (1.0 - gps[gp].x) * (1.0 + gps[gp].y) * (1.0 + gps[gp].z) * 0.125;
      sh[6][gp] = (1.0 + gps[gp].x) * (1.0 - gps[gp].y) * (1.0 + gps[gp].z) * 0.125;
      sh[7][gp] = (1.0 + gps[gp].x) * (1.0 + gps[gp].y) * (1.0 + gps[gp].z) * 0.125;
    }
    return sh;
  }

  constexpr Derivatives makeDerivatives() {
    Derivatives ds{};
    for (size_t gp = 0; gp < 8; gp++) {
      ds[0][0][gp] = -1.0 * (1.0 - gps[gp].y) * (1.0 - gps[gp].z) * 0.125;
      ds[0][1][gp] = (1.0 - gps[gp].x) * -1.0 * (1.0 - gps[gp].z) * 0.125;
      ds[0][2][gp] = (1.0 - gps[gp].x) * (1.0 - gps[gp].y) * -1.0 * 0.125;

      ds[1][0][gp] = -1.0 * (1.0 + gps[gp].y) * (1.0 - gps[gp].z) * 0.125;
      ds[1][1][gp] = (1.0 - gps[gp].x) * +1.0 * (1.0 - gps[gp].z) * 0.125;
      ds[1][2][gp] = (1.0 - gps[gp].x) * (1.0 + gps[gp].y) * -1.0 * 0.125;

      ds[2][0][gp] = +1.0 * (1.0 - gps[gp].y) * (1.0 - gps[gp].z) * 0.125;
      ds[2][1][gp] = (1.0 + gps[gp].x) * -1.0 * (1.0 - gps[gp].z) * 0.125;
      ds[2][2][gp] = (1.0 + gps[gp].x) * (1.0 - gps[gp].y) * -1.0 * 0.125;

      ds[3][0][gp] = +1.0 * (1.0 + gps[gp].y) * (1.0 - gps[gp].z) * 0.125;
      ds[3][1][gp] = (1.0 + gps[gp].x) * +1.0 * (1.0 - gps[gp].z) * 0.125;
      ds[3][2][gp] = (1.0 + gps[gp].x) * (1.0 + gps[gp].y) * -1.0 * 0.125;

      ds[4][0][gp] = -1.0 * (1.0 - gps[gp].y) * (1.0 + gps[gp].z) * 0.125;
      ds[4][1][gp] = (1.0 - gps[gp].x) * -1.0 * (1.0 + gps[gp].z) * 0.125;
      ds[4][2][gp] = (1.0 - gps[gp].x) * (1.0 - gps[gp].y) * +1.0 * 0.125;

      ds[5][0][gp] = -1.0 * (1.0 + gps[gp].y) * (1.0 + gps[gp].z) * 0.125;
      ds[5][1][gp] = (1.0 - gps[gp].x) * +1.0 * (1.0 + gps[gp].z) * 0.125;
      ds[5][2][gp] = (1.0 - gps[gp].x) * (1.0 + gps[gp].y) * +1.0 * 0.125;

      ds[6][0][gp] = +1.0 * (1.0 - gps[gp].y) * (1.0 + gps[gp].z) * 0.125;
      ds[6][1][gp] = (1.0 + gps[gp].x) * -1.0 * (1.0 + gps[gp].z) * 0.125;
      ds[6][2][gp] = (1.0 + gps[gp].x) * (1.0 - gps[gp].y) * +1.0 * 0.125;

      ds[7][0][gp] = +1.0 * (1.0 + gps[gp].y) * (1.0 + gps[gp].z) * 0.125;
      ds[7][1][gp] = (1.0 + gps[gp].x) * +1.0 * (1.0 + gps[gp].z) * 0.125;
      ds[7][2][gp] = (1.0 + gps[gp].x) * (1.0 + gps[gp].y) * +1.0 * 0.125;
    }
    return ds;
  }
};

#endif
