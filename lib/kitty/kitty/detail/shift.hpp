/* kitty: C++ truth table library
 * Copyright (C) 2017-2022  EPFL
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

/*!
  \file shift.hpp
  \brief Shift helper functions

  \author Mathias Soeken
*/

/*! \cond PRIVATE */
#pragma once

#include <cassert>
#include <cstdint>

namespace kitty
{
namespace detail
{

// This code has been auto-generated by Nikolaj Bjorner
// https://github.com/NikolajBjorner/z3/blob/17a256eca1a375f023bb98c3e4e728845d92dfb2/src/ast/expr_network.cpp#L22
inline uint64_t compute_shift( uint64_t x, unsigned code )
{
  switch ( code )
  {
#define v_x0 ( x & 1 )
#define v_x1 v_x0
  case 1:
    return v_x1;
#define v_x2 ( v_x1 | ( v_x1 << 1 ) )
  case 2:
    return v_x2;
#define v_x3 ( x & 3 )
#define v_x4 v_x3
  case 3:
    return v_x4;
#define v_x5 ( v_x2 | ( v_x2 << 2 ) )
  case 4:
    return v_x5;
#define v_x6 ( v_x4 | ( v_x4 << 2 ) )
  case 5:
    return v_x6;
#define v_x7 ( x & 2 )
#define v_x8 ( v_x7 << 1 )
#define v_x9 ( v_x8 | ( v_x8 << 1 ) )
#define v_x10 ( v_x2 | v_x9 )
  case 6:
    return v_x10;
#define v_x11 ( x & 15 )
#define v_x12 v_x11
  case 7:
    return v_x12;
#define v_x13 ( v_x5 | ( v_x5 << 4 ) )
  case 8:
    return v_x13;
#define v_x14 ( v_x6 | ( v_x6 << 4 ) )
  case 9:
    return v_x14;
#define v_x15 ( v_x10 | ( v_x10 << 4 ) )
  case 10:
    return v_x15;
#define v_x16 ( v_x12 | ( v_x12 << 4 ) )
  case 11:
    return v_x16;
#define v_x17 ( v_x7 << 3 )
#define v_x18 ( v_x17 | ( v_x17 << 1 ) )
#define v_x19 ( v_x18 | ( v_x18 << 2 ) )
#define v_x20 ( v_x5 | v_x19 )
  case 12:
    return v_x20;
#define v_x21 ( x & 12 )
#define v_x22 ( v_x21 << 2 )
#define v_x23 ( v_x22 | ( v_x22 << 2 ) )
#define v_x24 ( v_x6 | v_x23 )
  case 13:
    return v_x24;
#define v_x25 ( x & 4 )
#define v_x26 ( v_x25 << 2 )
#define v_x27 ( v_x26 | ( v_x26 << 1 ) )
#define v_x28 ( x & 8 )
#define v_x29 ( v_x28 << 3 )
#define v_x30 ( v_x29 | ( v_x29 << 1 ) )
#define v_x31 ( v_x27 | v_x30 )
#define v_x32 ( v_x10 | v_x31 )
  case 14:
    return v_x32;
#define v_x33 ( x & 255 )
#define v_x34 v_x33
  case 15:
    return v_x34;
#define v_x35 ( v_x13 | ( v_x13 << 8 ) )
  case 16:
    return v_x35;
#define v_x36 ( v_x14 | ( v_x14 << 8 ) )
  case 17:
    return v_x36;
#define v_x37 ( v_x15 | ( v_x15 << 8 ) )
  case 18:
    return v_x37;
#define v_x38 ( v_x16 | ( v_x16 << 8 ) )
  case 19:
    return v_x38;
#define v_x39 ( v_x20 | ( v_x20 << 8 ) )
  case 20:
    return v_x39;
#define v_x40 ( v_x24 | ( v_x24 << 8 ) )
  case 21:
    return v_x40;
#define v_x41 ( v_x32 | ( v_x32 << 8 ) )
  case 22:
    return v_x41;
#define v_x42 ( v_x34 | ( v_x34 << 8 ) )
  case 23:
    return v_x42;
#define v_x43 ( v_x7 << 7 )
#define v_x44 ( v_x43 | ( v_x43 << 1 ) )
#define v_x45 ( v_x44 | ( v_x44 << 2 ) )
#define v_x46 ( v_x45 | ( v_x45 << 4 ) )
#define v_x47 ( v_x13 | v_x46 )
  case 24:
    return v_x47;
#define v_x48 ( v_x21 << 6 )
#define v_x49 ( v_x48 | ( v_x48 << 2 ) )
#define v_x50 ( v_x49 | ( v_x49 << 4 ) )
#define v_x51 ( v_x14 | v_x50 )
  case 25:
    return v_x51;
#define v_x52 ( v_x25 << 6 )
#define v_x53 ( v_x52 | ( v_x52 << 1 ) )
#define v_x54 ( v_x28 << 7 )
#define v_x55 ( v_x54 | ( v_x54 << 1 ) )
#define v_x56 ( v_x53 | v_x55 )
#define v_x57 ( v_x56 | ( v_x56 << 4 ) )
#define v_x58 ( v_x15 | v_x57 )
  case 26:
    return v_x58;
#define v_x59 ( x & 240 )
#define v_x60 ( v_x59 << 4 )
#define v_x61 ( v_x60 | ( v_x60 << 4 ) )
#define v_x62 ( v_x16 | v_x61 )
  case 27:
    return v_x62;
#define v_x63 ( v_x53 | ( v_x53 << 2 ) )
#define v_x64 ( v_x28 << 9 )
#define v_x65 ( v_x64 | ( v_x64 << 1 ) )
#define v_x66 ( v_x65 | ( v_x65 << 2 ) )
#define v_x67 ( v_x63 | v_x66 )
#define v_x68 ( v_x20 | v_x67 )
  case 28:
    return v_x68;
#define v_x69 ( x & 48 )
#define v_x70 ( v_x69 << 4 )
#define v_x71 ( v_x70 | ( v_x70 << 2 ) )
#define v_x72 ( x & 192 )
#define v_x73 ( v_x72 << 6 )
#define v_x74 ( v_x73 | ( v_x73 << 2 ) )
#define v_x75 ( v_x71 | v_x74 )
#define v_x76 ( v_x24 | v_x75 )
  case 29:
    return v_x76;
#define v_x77 ( x & 16 )
#define v_x78 ( v_x77 << 4 )
#define v_x79 ( v_x78 | ( v_x78 << 1 ) )
#define v_x80 ( x & 32 )
#define v_x81 ( v_x80 << 5 )
#define v_x82 ( v_x81 | ( v_x81 << 1 ) )
#define v_x83 ( v_x79 | v_x82 )
#define v_x84 ( x & 64 )
#define v_x85 ( v_x84 << 6 )
#define v_x86 ( v_x85 | ( v_x85 << 1 ) )
#define v_x87 ( x & 128 )
#define v_x88 ( v_x87 << 7 )
#define v_x89 ( v_x88 | ( v_x88 << 1 ) )
#define v_x90 ( v_x86 | v_x89 )
#define v_x91 ( v_x83 | v_x90 )
#define v_x92 ( v_x32 | v_x91 )
  case 30:
    return v_x92;
#define v_x93 ( x & 65535 )
#define v_x94 v_x93
  case 31:
    return v_x94;
#define v_x95 ( v_x35 | ( v_x35 << 16 ) )
  case 32:
    return v_x95;
#define v_x96 ( v_x36 | ( v_x36 << 16 ) )
  case 33:
    return v_x96;
#define v_x97 ( v_x37 | ( v_x37 << 16 ) )
  case 34:
    return v_x97;
#define v_x98 ( v_x38 | ( v_x38 << 16 ) )
  case 35:
    return v_x98;
#define v_x99 ( v_x39 | ( v_x39 << 16 ) )
  case 36:
    return v_x99;
#define v_x100 ( v_x40 | ( v_x40 << 16 ) )
  case 37:
    return v_x100;
#define v_x101 ( v_x41 | ( v_x41 << 16 ) )
  case 38:
    return v_x101;
#define v_x102 ( v_x42 | ( v_x42 << 16 ) )
  case 39:
    return v_x102;
#define v_x103 ( v_x47 | ( v_x47 << 16 ) )
  case 40:
    return v_x103;
#define v_x104 ( v_x51 | ( v_x51 << 16 ) )
  case 41:
    return v_x104;
#define v_x105 ( v_x58 | ( v_x58 << 16 ) )
  case 42:
    return v_x105;
#define v_x106 ( v_x62 | ( v_x62 << 16 ) )
  case 43:
    return v_x106;
#define v_x107 ( v_x68 | ( v_x68 << 16 ) )
  case 44:
    return v_x107;
#define v_x108 ( v_x76 | ( v_x76 << 16 ) )
  case 45:
    return v_x108;
#define v_x109 ( v_x92 | ( v_x92 << 16 ) )
  case 46:
    return v_x109;
#define v_x110 ( v_x94 | ( v_x94 << 16 ) )
  case 47:
    return v_x110;
#define v_x111 ( v_x7 << 15 )
#define v_x112 ( v_x111 | ( v_x111 << 1 ) )
#define v_x113 ( v_x112 | ( v_x112 << 2 ) )
#define v_x114 ( v_x113 | ( v_x113 << 4 ) )
#define v_x115 ( v_x114 | ( v_x114 << 8 ) )
#define v_x116 ( v_x35 | v_x115 )
  case 48:
    return v_x116;
#define v_x117 ( v_x21 << 14 )
#define v_x118 ( v_x117 | ( v_x117 << 2 ) )
#define v_x119 ( v_x118 | ( v_x118 << 4 ) )
#define v_x120 ( v_x119 | ( v_x119 << 8 ) )
#define v_x121 ( v_x36 | v_x120 )
  case 49:
    return v_x121;
#define v_x122 ( v_x25 << 14 )
#define v_x123 ( v_x122 | ( v_x122 << 1 ) )
#define v_x124 ( v_x28 << 15 )
#define v_x125 ( v_x124 | ( v_x124 << 1 ) )
#define v_x126 ( v_x123 | v_x125 )
#define v_x127 ( v_x126 | ( v_x126 << 4 ) )
#define v_x128 ( v_x127 | ( v_x127 << 8 ) )
#define v_x129 ( v_x37 | v_x128 )
  case 50:
    return v_x129;
#define v_x130 ( v_x59 << 12 )
#define v_x131 ( v_x130 | ( v_x130 << 4 ) )
#define v_x132 ( v_x131 | ( v_x131 << 8 ) )
#define v_x133 ( v_x38 | v_x132 )
  case 51:
    return v_x133;
#define v_x134 ( v_x123 | ( v_x123 << 2 ) )
#define v_x135 ( v_x28 << 17 )
#define v_x136 ( v_x135 | ( v_x135 << 1 ) )
#define v_x137 ( v_x136 | ( v_x136 << 2 ) )
#define v_x138 ( v_x134 | v_x137 )
#define v_x139 ( v_x138 | ( v_x138 << 8 ) )
#define v_x140 ( v_x39 | v_x139 )
  case 52:
    return v_x140;
#define v_x141 ( v_x69 << 12 )
#define v_x142 ( v_x141 | ( v_x141 << 2 ) )
#define v_x143 ( v_x72 << 14 )
#define v_x144 ( v_x143 | ( v_x143 << 2 ) )
#define v_x145 ( v_x142 | v_x144 )
#define v_x146 ( v_x145 | ( v_x145 << 8 ) )
#define v_x147 ( v_x40 | v_x146 )
  case 53:
    return v_x147;
#define v_x148 ( v_x77 << 12 )
#define v_x149 ( v_x148 | ( v_x148 << 1 ) )
#define v_x150 ( v_x80 << 13 )
#define v_x151 ( v_x150 | ( v_x150 << 1 ) )
#define v_x152 ( v_x149 | v_x151 )
#define v_x153 ( v_x84 << 14 )
#define v_x154 ( v_x153 | ( v_x153 << 1 ) )
#define v_x155 ( v_x87 << 15 )
#define v_x156 ( v_x155 | ( v_x155 << 1 ) )
#define v_x157 ( v_x154 | v_x156 )
#define v_x158 ( v_x152 | v_x157 )
#define v_x159 ( v_x158 | ( v_x158 << 8 ) )
#define v_x160 ( v_x41 | v_x159 )
  case 54:
    return v_x160;
#define v_x161 ( x & 65280 )
#define v_x162 ( v_x161 << 8 )
#define v_x163 ( v_x162 | ( v_x162 << 8 ) )
#define v_x164 ( v_x42 | v_x163 )
  case 55:
    return v_x164;
#define v_x165 ( v_x134 | ( v_x134 << 4 ) )
#define v_x166 ( v_x28 << 21 )
#define v_x167 ( v_x166 | ( v_x166 << 1 ) )
#define v_x168 ( v_x167 | ( v_x167 << 2 ) )
#define v_x169 ( v_x168 | ( v_x168 << 4 ) )
#define v_x170 ( v_x165 | v_x169 )
#define v_x171 ( v_x47 | v_x170 )
  case 56:
    return v_x171;
#define v_x172 ( v_x142 | ( v_x142 << 4 ) )
#define v_x173 ( v_x72 << 18 )
#define v_x174 ( v_x173 | ( v_x173 << 2 ) )
#define v_x175 ( v_x174 | ( v_x174 << 4 ) )
#define v_x176 ( v_x172 | v_x175 )
#define v_x177 ( v_x51 | v_x176 )
  case 57:
    return v_x177;
#define v_x178 ( v_x152 | ( v_x152 << 4 ) )
#define v_x179 ( v_x84 << 18 )
#define v_x180 ( v_x179 | ( v_x179 << 1 ) )
#define v_x181 ( v_x87 << 19 )
#define v_x182 ( v_x181 | ( v_x181 << 1 ) )
#define v_x183 ( v_x180 | v_x182 )
#define v_x184 ( v_x183 | ( v_x183 << 4 ) )
#define v_x185 ( v_x178 | v_x184 )
#define v_x186 ( v_x58 | v_x185 )
  case 58:
    return v_x186;
#define v_x187 ( x & 3840 )
#define v_x188 ( v_x187 << 8 )
#define v_x189 ( v_x188 | ( v_x188 << 4 ) )
#define v_x190 ( x & 61440 )
#define v_x191 ( v_x190 << 12 )
#define v_x192 ( v_x191 | ( v_x191 << 4 ) )
#define v_x193 ( v_x189 | v_x192 )
#define v_x194 ( v_x62 | v_x193 )
  case 59:
    return v_x194;
#define v_x195 ( v_x149 | ( v_x149 << 2 ) )
#define v_x196 ( v_x80 << 15 )
#define v_x197 ( v_x196 | ( v_x196 << 1 ) )
#define v_x198 ( v_x197 | ( v_x197 << 2 ) )
#define v_x199 ( v_x195 | v_x198 )
#define v_x200 ( v_x180 | ( v_x180 << 2 ) )
#define v_x201 ( v_x87 << 21 )
#define v_x202 ( v_x201 | ( v_x201 << 1 ) )
#define v_x203 ( v_x202 | ( v_x202 << 2 ) )
#define v_x204 ( v_x200 | v_x203 )
#define v_x205 ( v_x199 | v_x204 )
#define v_x206 ( v_x68 | v_x205 )
  case 60:
    return v_x206;
#define v_x207 ( x & 768 )
#define v_x208 ( v_x207 << 8 )
#define v_x209 ( v_x208 | ( v_x208 << 2 ) )
#define v_x210 ( x & 3072 )
#define v_x211 ( v_x210 << 10 )
#define v_x212 ( v_x211 | ( v_x211 << 2 ) )
#define v_x213 ( v_x209 | v_x212 )
#define v_x214 ( x & 12288 )
#define v_x215 ( v_x214 << 12 )
#define v_x216 ( v_x215 | ( v_x215 << 2 ) )
#define v_x217 ( x & 49152 )
#define v_x218 ( v_x217 << 14 )
#define v_x219 ( v_x218 | ( v_x218 << 2 ) )
#define v_x220 ( v_x216 | v_x219 )
#define v_x221 ( v_x213 | v_x220 )
#define v_x222 ( v_x76 | v_x221 )
  case 61:
    return v_x222;
#define v_x223 ( x & 256 )
#define v_x224 ( v_x223 << 8 )
#define v_x225 ( v_x224 | ( v_x224 << 1 ) )
#define v_x226 ( x & 512 )
#define v_x227 ( v_x226 << 9 )
#define v_x228 ( v_x227 | ( v_x227 << 1 ) )
#define v_x229 ( v_x225 | v_x228 )
#define v_x230 ( x & 1024 )
#define v_x231 ( v_x230 << 10 )
#define v_x232 ( v_x231 | ( v_x231 << 1 ) )
#define v_x233 ( x & 2048 )
#define v_x234 ( v_x233 << 11 )
#define v_x235 ( v_x234 | ( v_x234 << 1 ) )
#define v_x236 ( v_x232 | v_x235 )
#define v_x237 ( v_x229 | v_x236 )
#define v_x238 ( x & 4096 )
#define v_x239 ( v_x238 << 12 )
#define v_x240 ( v_x239 | ( v_x239 << 1 ) )
#define v_x241 ( x & 8192 )
#define v_x242 ( v_x241 << 13 )
#define v_x243 ( v_x242 | ( v_x242 << 1 ) )
#define v_x244 ( v_x240 | v_x243 )
#define v_x245 ( x & 16384 )
#define v_x246 ( v_x245 << 14 )
#define v_x247 ( v_x246 | ( v_x246 << 1 ) )
#define v_x248 ( x & 32768 )
#define v_x249 ( v_x248 << 15 )
#define v_x250 ( v_x249 | ( v_x249 << 1 ) )
#define v_x251 ( v_x247 | v_x250 )
#define v_x252 ( v_x244 | v_x251 )
#define v_x253 ( v_x237 | v_x252 )
#define v_x254 ( v_x92 | v_x253 )
  case 62:
    return v_x254;
#define v_x255 ( x & 4294967295 )
#define v_x256 v_x255
  case 63:
    return v_x256;
#define v_x257 ( v_x95 | ( v_x95 << 32 ) )
  case 64:
    return v_x257;
#define v_x258 ( v_x96 | ( v_x96 << 32 ) )
  case 65:
    return v_x258;
#define v_x259 ( v_x97 | ( v_x97 << 32 ) )
  case 66:
    return v_x259;
#define v_x260 ( v_x98 | ( v_x98 << 32 ) )
  case 67:
    return v_x260;
#define v_x261 ( v_x99 | ( v_x99 << 32 ) )
  case 68:
    return v_x261;
#define v_x262 ( v_x100 | ( v_x100 << 32 ) )
  case 69:
    return v_x262;
#define v_x263 ( v_x101 | ( v_x101 << 32 ) )
  case 70:
    return v_x263;
#define v_x264 ( v_x102 | ( v_x102 << 32 ) )
  case 71:
    return v_x264;
#define v_x265 ( v_x103 | ( v_x103 << 32 ) )
  case 72:
    return v_x265;
#define v_x266 ( v_x104 | ( v_x104 << 32 ) )
  case 73:
    return v_x266;
#define v_x267 ( v_x105 | ( v_x105 << 32 ) )
  case 74:
    return v_x267;
#define v_x268 ( v_x106 | ( v_x106 << 32 ) )
  case 75:
    return v_x268;
#define v_x269 ( v_x107 | ( v_x107 << 32 ) )
  case 76:
    return v_x269;
#define v_x270 ( v_x108 | ( v_x108 << 32 ) )
  case 77:
    return v_x270;
#define v_x271 ( v_x109 | ( v_x109 << 32 ) )
  case 78:
    return v_x271;
#define v_x272 ( v_x110 | ( v_x110 << 32 ) )
  case 79:
    return v_x272;
#define v_x273 ( v_x116 | ( v_x116 << 32 ) )
  case 80:
    return v_x273;
#define v_x274 ( v_x121 | ( v_x121 << 32 ) )
  case 81:
    return v_x274;
#define v_x275 ( v_x129 | ( v_x129 << 32 ) )
  case 82:
    return v_x275;
#define v_x276 ( v_x133 | ( v_x133 << 32 ) )
  case 83:
    return v_x276;
#define v_x277 ( v_x140 | ( v_x140 << 32 ) )
  case 84:
    return v_x277;
#define v_x278 ( v_x147 | ( v_x147 << 32 ) )
  case 85:
    return v_x278;
#define v_x279 ( v_x160 | ( v_x160 << 32 ) )
  case 86:
    return v_x279;
#define v_x280 ( v_x164 | ( v_x164 << 32 ) )
  case 87:
    return v_x280;
#define v_x281 ( v_x171 | ( v_x171 << 32 ) )
  case 88:
    return v_x281;
#define v_x282 ( v_x177 | ( v_x177 << 32 ) )
  case 89:
    return v_x282;
#define v_x283 ( v_x186 | ( v_x186 << 32 ) )
  case 90:
    return v_x283;
#define v_x284 ( v_x194 | ( v_x194 << 32 ) )
  case 91:
    return v_x284;
#define v_x285 ( v_x206 | ( v_x206 << 32 ) )
  case 92:
    return v_x285;
#define v_x286 ( v_x222 | ( v_x222 << 32 ) )
  case 93:
    return v_x286;
#define v_x287 ( v_x254 | ( v_x254 << 32 ) )
  case 94:
    return v_x287;
#define v_x288 ( v_x256 | ( v_x256 << 32 ) )
  case 95:
    return v_x288;
#define v_x289 ( v_x7 << 31 )
#define v_x290 ( v_x289 | ( v_x289 << 1 ) )
#define v_x291 ( v_x290 | ( v_x290 << 2 ) )
#define v_x292 ( v_x291 | ( v_x291 << 4 ) )
#define v_x293 ( v_x292 | ( v_x292 << 8 ) )
#define v_x294 ( v_x293 | ( v_x293 << 16 ) )
#define v_x295 ( v_x95 | v_x294 )
  case 96:
    return v_x295;
#define v_x296 ( v_x21 << 30 )
#define v_x297 ( v_x296 | ( v_x296 << 2 ) )
#define v_x298 ( v_x297 | ( v_x297 << 4 ) )
#define v_x299 ( v_x298 | ( v_x298 << 8 ) )
#define v_x300 ( v_x299 | ( v_x299 << 16 ) )
#define v_x301 ( v_x96 | v_x300 )
  case 97:
    return v_x301;
#define v_x302 ( v_x25 << 30 )
#define v_x303 ( v_x302 | ( v_x302 << 1 ) )
#define v_x304 ( v_x28 << 31 )
#define v_x305 ( v_x304 | ( v_x304 << 1 ) )
#define v_x306 ( v_x303 | v_x305 )
#define v_x307 ( v_x306 | ( v_x306 << 4 ) )
#define v_x308 ( v_x307 | ( v_x307 << 8 ) )
#define v_x309 ( v_x308 | ( v_x308 << 16 ) )
#define v_x310 ( v_x97 | v_x309 )
  case 98:
    return v_x310;
#define v_x311 ( v_x59 << 28 )
#define v_x312 ( v_x311 | ( v_x311 << 4 ) )
#define v_x313 ( v_x312 | ( v_x312 << 8 ) )
#define v_x314 ( v_x313 | ( v_x313 << 16 ) )
#define v_x315 ( v_x98 | v_x314 )
  case 99:
    return v_x315;
#define v_x316 ( v_x303 | ( v_x303 << 2 ) )
#define v_x317 ( v_x28 << 33 )
#define v_x318 ( v_x317 | ( v_x317 << 1 ) )
#define v_x319 ( v_x318 | ( v_x318 << 2 ) )
#define v_x320 ( v_x316 | v_x319 )
#define v_x321 ( v_x320 | ( v_x320 << 8 ) )
#define v_x322 ( v_x321 | ( v_x321 << 16 ) )
#define v_x323 ( v_x99 | v_x322 )
  case 100:
    return v_x323;
#define v_x324 ( v_x69 << 28 )
#define v_x325 ( v_x324 | ( v_x324 << 2 ) )
#define v_x326 ( v_x72 << 30 )
#define v_x327 ( v_x326 | ( v_x326 << 2 ) )
#define v_x328 ( v_x325 | v_x327 )
#define v_x329 ( v_x328 | ( v_x328 << 8 ) )
#define v_x330 ( v_x329 | ( v_x329 << 16 ) )
#define v_x331 ( v_x100 | v_x330 )
  case 101:
    return v_x331;
#define v_x332 ( v_x77 << 28 )
#define v_x333 ( v_x332 | ( v_x332 << 1 ) )
#define v_x334 ( v_x80 << 29 )
#define v_x335 ( v_x334 | ( v_x334 << 1 ) )
#define v_x336 ( v_x333 | v_x335 )
#define v_x337 ( v_x84 << 30 )
#define v_x338 ( v_x337 | ( v_x337 << 1 ) )
#define v_x339 ( v_x87 << 31 )
#define v_x340 ( v_x339 | ( v_x339 << 1 ) )
#define v_x341 ( v_x338 | v_x340 )
#define v_x342 ( v_x336 | v_x341 )
#define v_x343 ( v_x342 | ( v_x342 << 8 ) )
#define v_x344 ( v_x343 | ( v_x343 << 16 ) )
#define v_x345 ( v_x101 | v_x344 )
  case 102:
    return v_x345;
#define v_x346 ( v_x161 << 24 )
#define v_x347 ( v_x346 | ( v_x346 << 8 ) )
#define v_x348 ( v_x347 | ( v_x347 << 16 ) )
#define v_x349 ( v_x102 | v_x348 )
  case 103:
    return v_x349;
#define v_x350 ( v_x316 | ( v_x316 << 4 ) )
#define v_x351 ( v_x28 << 37 )
#define v_x352 ( v_x351 | ( v_x351 << 1 ) )
#define v_x353 ( v_x352 | ( v_x352 << 2 ) )
#define v_x354 ( v_x353 | ( v_x353 << 4 ) )
#define v_x355 ( v_x350 | v_x354 )
#define v_x356 ( v_x355 | ( v_x355 << 16 ) )
#define v_x357 ( v_x103 | v_x356 )
  case 104:
    return v_x357;
#define v_x358 ( v_x325 | ( v_x325 << 4 ) )
#define v_x359 ( v_x72 << 34 )
#define v_x360 ( v_x359 | ( v_x359 << 2 ) )
#define v_x361 ( v_x360 | ( v_x360 << 4 ) )
#define v_x362 ( v_x358 | v_x361 )
#define v_x363 ( v_x362 | ( v_x362 << 16 ) )
#define v_x364 ( v_x104 | v_x363 )
  case 105:
    return v_x364;
#define v_x365 ( v_x336 | ( v_x336 << 4 ) )
#define v_x366 ( v_x84 << 34 )
#define v_x367 ( v_x366 | ( v_x366 << 1 ) )
#define v_x368 ( v_x87 << 35 )
#define v_x369 ( v_x368 | ( v_x368 << 1 ) )
#define v_x370 ( v_x367 | v_x369 )
#define v_x371 ( v_x370 | ( v_x370 << 4 ) )
#define v_x372 ( v_x365 | v_x371 )
#define v_x373 ( v_x372 | ( v_x372 << 16 ) )
#define v_x374 ( v_x105 | v_x373 )
  case 106:
    return v_x374;
#define v_x375 ( v_x187 << 24 )
#define v_x376 ( v_x375 | ( v_x375 << 4 ) )
#define v_x377 ( v_x190 << 28 )
#define v_x378 ( v_x377 | ( v_x377 << 4 ) )
#define v_x379 ( v_x376 | v_x378 )
#define v_x380 ( v_x379 | ( v_x379 << 16 ) )
#define v_x381 ( v_x106 | v_x380 )
  case 107:
    return v_x381;
#define v_x382 ( v_x333 | ( v_x333 << 2 ) )
#define v_x383 ( v_x80 << 31 )
#define v_x384 ( v_x383 | ( v_x383 << 1 ) )
#define v_x385 ( v_x384 | ( v_x384 << 2 ) )
#define v_x386 ( v_x382 | v_x385 )
#define v_x387 ( v_x367 | ( v_x367 << 2 ) )
#define v_x388 ( v_x87 << 37 )
#define v_x389 ( v_x388 | ( v_x388 << 1 ) )
#define v_x390 ( v_x389 | ( v_x389 << 2 ) )
#define v_x391 ( v_x387 | v_x390 )
#define v_x392 ( v_x386 | v_x391 )
#define v_x393 ( v_x392 | ( v_x392 << 16 ) )
#define v_x394 ( v_x107 | v_x393 )
  case 108:
    return v_x394;
#define v_x395 ( v_x207 << 24 )
#define v_x396 ( v_x395 | ( v_x395 << 2 ) )
#define v_x397 ( v_x210 << 26 )
#define v_x398 ( v_x397 | ( v_x397 << 2 ) )
#define v_x399 ( v_x396 | v_x398 )
#define v_x400 ( v_x214 << 28 )
#define v_x401 ( v_x400 | ( v_x400 << 2 ) )
#define v_x402 ( v_x217 << 30 )
#define v_x403 ( v_x402 | ( v_x402 << 2 ) )
#define v_x404 ( v_x401 | v_x403 )
#define v_x405 ( v_x399 | v_x404 )
#define v_x406 ( v_x405 | ( v_x405 << 16 ) )
#define v_x407 ( v_x108 | v_x406 )
  case 109:
    return v_x407;
#define v_x408 ( v_x223 << 24 )
#define v_x409 ( v_x408 | ( v_x408 << 1 ) )
#define v_x410 ( v_x226 << 25 )
#define v_x411 ( v_x410 | ( v_x410 << 1 ) )
#define v_x412 ( v_x409 | v_x411 )
#define v_x413 ( v_x230 << 26 )
#define v_x414 ( v_x413 | ( v_x413 << 1 ) )
#define v_x415 ( v_x233 << 27 )
#define v_x416 ( v_x415 | ( v_x415 << 1 ) )
#define v_x417 ( v_x414 | v_x416 )
#define v_x418 ( v_x412 | v_x417 )
#define v_x419 ( v_x238 << 28 )
#define v_x420 ( v_x419 | ( v_x419 << 1 ) )
#define v_x421 ( v_x241 << 29 )
#define v_x422 ( v_x421 | ( v_x421 << 1 ) )
#define v_x423 ( v_x420 | v_x422 )
#define v_x424 ( v_x245 << 30 )
#define v_x425 ( v_x424 | ( v_x424 << 1 ) )
#define v_x426 ( v_x248 << 31 )
#define v_x427 ( v_x426 | ( v_x426 << 1 ) )
#define v_x428 ( v_x425 | v_x427 )
#define v_x429 ( v_x423 | v_x428 )
#define v_x430 ( v_x418 | v_x429 )
#define v_x431 ( v_x430 | ( v_x430 << 16 ) )
#define v_x432 ( v_x109 | v_x431 )
  case 110:
    return v_x432;
#define v_x433 ( x & 4294901760 )
#define v_x434 ( v_x433 << 16 )
#define v_x435 ( v_x434 | ( v_x434 << 16 ) )
#define v_x436 ( v_x110 | v_x435 )
  case 111:
    return v_x436;
#define v_x437 ( v_x350 | ( v_x350 << 8 ) )
#define v_x438 ( v_x28 << 45 )
#define v_x439 ( v_x438 | ( v_x438 << 1 ) )
#define v_x440 ( v_x439 | ( v_x439 << 2 ) )
#define v_x441 ( v_x440 | ( v_x440 << 4 ) )
#define v_x442 ( v_x441 | ( v_x441 << 8 ) )
#define v_x443 ( v_x437 | v_x442 )
#define v_x444 ( v_x116 | v_x443 )
  case 112:
    return v_x444;
#define v_x445 ( v_x358 | ( v_x358 << 8 ) )
#define v_x446 ( v_x72 << 42 )
#define v_x447 ( v_x446 | ( v_x446 << 2 ) )
#define v_x448 ( v_x447 | ( v_x447 << 4 ) )
#define v_x449 ( v_x448 | ( v_x448 << 8 ) )
#define v_x450 ( v_x445 | v_x449 )
#define v_x451 ( v_x121 | v_x450 )
  case 113:
    return v_x451;
#define v_x452 ( v_x365 | ( v_x365 << 8 ) )
#define v_x453 ( v_x84 << 42 )
#define v_x454 ( v_x453 | ( v_x453 << 1 ) )
#define v_x455 ( v_x87 << 43 )
#define v_x456 ( v_x455 | ( v_x455 << 1 ) )
#define v_x457 ( v_x454 | v_x456 )
#define v_x458 ( v_x457 | ( v_x457 << 4 ) )
#define v_x459 ( v_x458 | ( v_x458 << 8 ) )
#define v_x460 ( v_x452 | v_x459 )
#define v_x461 ( v_x129 | v_x460 )
  case 114:
    return v_x461;
#define v_x462 ( v_x376 | ( v_x376 << 8 ) )
#define v_x463 ( v_x190 << 36 )
#define v_x464 ( v_x463 | ( v_x463 << 4 ) )
#define v_x465 ( v_x464 | ( v_x464 << 8 ) )
#define v_x466 ( v_x462 | v_x465 )
#define v_x467 ( v_x133 | v_x466 )
  case 115:
    return v_x467;
#define v_x468 ( v_x386 | ( v_x386 << 8 ) )
#define v_x469 ( v_x454 | ( v_x454 << 2 ) )
#define v_x470 ( v_x87 << 45 )
#define v_x471 ( v_x470 | ( v_x470 << 1 ) )
#define v_x472 ( v_x471 | ( v_x471 << 2 ) )
#define v_x473 ( v_x469 | v_x472 )
#define v_x474 ( v_x473 | ( v_x473 << 8 ) )
#define v_x475 ( v_x468 | v_x474 )
#define v_x476 ( v_x140 | v_x475 )
  case 116:
    return v_x476;
#define v_x477 ( v_x399 | ( v_x399 << 8 ) )
#define v_x478 ( v_x214 << 36 )
#define v_x479 ( v_x478 | ( v_x478 << 2 ) )
#define v_x480 ( v_x217 << 38 )
#define v_x481 ( v_x480 | ( v_x480 << 2 ) )
#define v_x482 ( v_x479 | v_x481 )
#define v_x483 ( v_x482 | ( v_x482 << 8 ) )
#define v_x484 ( v_x477 | v_x483 )
#define v_x485 ( v_x147 | v_x484 )
  case 117:
    return v_x485;
#define v_x486 ( v_x418 | ( v_x418 << 8 ) )
#define v_x487 ( v_x238 << 36 )
#define v_x488 ( v_x487 | ( v_x487 << 1 ) )
#define v_x489 ( v_x241 << 37 )
#define v_x490 ( v_x489 | ( v_x489 << 1 ) )
#define v_x491 ( v_x488 | v_x490 )
#define v_x492 ( v_x245 << 38 )
#define v_x493 ( v_x492 | ( v_x492 << 1 ) )
#define v_x494 ( v_x248 << 39 )
#define v_x495 ( v_x494 | ( v_x494 << 1 ) )
#define v_x496 ( v_x493 | v_x495 )
#define v_x497 ( v_x491 | v_x496 )
#define v_x498 ( v_x497 | ( v_x497 << 8 ) )
#define v_x499 ( v_x486 | v_x498 )
#define v_x500 ( v_x160 | v_x499 )
  case 118:
    return v_x500;
#define v_x501 ( x & 16711680 )
#define v_x502 ( v_x501 << 16 )
#define v_x503 ( v_x502 | ( v_x502 << 8 ) )
#define v_x504 ( x & 4278190080 )
#define v_x505 ( v_x504 << 24 )
#define v_x506 ( v_x505 | ( v_x505 << 8 ) )
#define v_x507 ( v_x503 | v_x506 )
#define v_x508 ( v_x164 | v_x507 )
  case 119:
    return v_x508;
#define v_x509 ( v_x382 | ( v_x382 << 4 ) )
#define v_x510 ( v_x80 << 35 )
#define v_x511 ( v_x510 | ( v_x510 << 1 ) )
#define v_x512 ( v_x511 | ( v_x511 << 2 ) )
#define v_x513 ( v_x512 | ( v_x512 << 4 ) )
#define v_x514 ( v_x509 | v_x513 )
#define v_x515 ( v_x469 | ( v_x469 << 4 ) )
#define v_x516 ( v_x87 << 49 )
#define v_x517 ( v_x516 | ( v_x516 << 1 ) )
#define v_x518 ( v_x517 | ( v_x517 << 2 ) )
#define v_x519 ( v_x518 | ( v_x518 << 4 ) )
#define v_x520 ( v_x515 | v_x519 )
#define v_x521 ( v_x514 | v_x520 )
#define v_x522 ( v_x171 | v_x521 )
  case 120:
    return v_x522;
#define v_x523 ( v_x396 | ( v_x396 << 4 ) )
#define v_x524 ( v_x210 << 30 )
#define v_x525 ( v_x524 | ( v_x524 << 2 ) )
#define v_x526 ( v_x525 | ( v_x525 << 4 ) )
#define v_x527 ( v_x523 | v_x526 )
#define v_x528 ( v_x479 | ( v_x479 << 4 ) )
#define v_x529 ( v_x217 << 42 )
#define v_x530 ( v_x529 | ( v_x529 << 2 ) )
#define v_x531 ( v_x530 | ( v_x530 << 4 ) )
#define v_x532 ( v_x528 | v_x531 )
#define v_x533 ( v_x527 | v_x532 )
#define v_x534 ( v_x177 | v_x533 )
  case 121:
    return v_x534;
#define v_x535 ( v_x412 | ( v_x412 << 4 ) )
#define v_x536 ( v_x230 << 30 )
#define v_x537 ( v_x536 | ( v_x536 << 1 ) )
#define v_x538 ( v_x233 << 31 )
#define v_x539 ( v_x538 | ( v_x538 << 1 ) )
#define v_x540 ( v_x537 | v_x539 )
#define v_x541 ( v_x540 | ( v_x540 << 4 ) )
#define v_x542 ( v_x535 | v_x541 )
#define v_x543 ( v_x491 | ( v_x491 << 4 ) )
#define v_x544 ( v_x245 << 42 )
#define v_x545 ( v_x544 | ( v_x544 << 1 ) )
#define v_x546 ( v_x248 << 43 )
#define v_x547 ( v_x546 | ( v_x546 << 1 ) )
#define v_x548 ( v_x545 | v_x547 )
#define v_x549 ( v_x548 | ( v_x548 << 4 ) )
#define v_x550 ( v_x543 | v_x549 )
#define v_x551 ( v_x542 | v_x550 )
#define v_x552 ( v_x186 | v_x551 )
  case 122:
    return v_x552;
#define v_x553 ( x & 983040 )
#define v_x554 ( v_x553 << 16 )
#define v_x555 ( v_x554 | ( v_x554 << 4 ) )
#define v_x556 ( x & 15728640 )
#define v_x557 ( v_x556 << 20 )
#define v_x558 ( v_x557 | ( v_x557 << 4 ) )
#define v_x559 ( v_x555 | v_x558 )
#define v_x560 ( x & 251658240 )
#define v_x561 ( v_x560 << 24 )
#define v_x562 ( v_x561 | ( v_x561 << 4 ) )
#define v_x563 ( x & 4026531840 )
#define v_x564 ( v_x563 << 28 )
#define v_x565 ( v_x564 | ( v_x564 << 4 ) )
#define v_x566 ( v_x562 | v_x565 )
#define v_x567 ( v_x559 | v_x566 )
#define v_x568 ( v_x194 | v_x567 )
  case 123:
    return v_x568;
#define v_x569 ( v_x409 | ( v_x409 << 2 ) )
#define v_x570 ( v_x226 << 27 )
#define v_x571 ( v_x570 | ( v_x570 << 1 ) )
#define v_x572 ( v_x571 | ( v_x571 << 2 ) )
#define v_x573 ( v_x569 | v_x572 )
#define v_x574 ( v_x537 | ( v_x537 << 2 ) )
#define v_x575 ( v_x233 << 33 )
#define v_x576 ( v_x575 | ( v_x575 << 1 ) )
#define v_x577 ( v_x576 | ( v_x576 << 2 ) )
#define v_x578 ( v_x574 | v_x577 )
#define v_x579 ( v_x573 | v_x578 )
#define v_x580 ( v_x488 | ( v_x488 << 2 ) )
#define v_x581 ( v_x241 << 39 )
#define v_x582 ( v_x581 | ( v_x581 << 1 ) )
#define v_x583 ( v_x582 | ( v_x582 << 2 ) )
#define v_x584 ( v_x580 | v_x583 )
#define v_x585 ( v_x545 | ( v_x545 << 2 ) )
#define v_x586 ( v_x248 << 45 )
#define v_x587 ( v_x586 | ( v_x586 << 1 ) )
#define v_x588 ( v_x587 | ( v_x587 << 2 ) )
#define v_x589 ( v_x585 | v_x588 )
#define v_x590 ( v_x584 | v_x589 )
#define v_x591 ( v_x579 | v_x590 )
#define v_x592 ( v_x206 | v_x591 )
  case 124:
    return v_x592;
#define v_x593 ( x & 196608 )
#define v_x594 ( v_x593 << 16 )
#define v_x595 ( v_x594 | ( v_x594 << 2 ) )
#define v_x596 ( x & 786432 )
#define v_x597 ( v_x596 << 18 )
#define v_x598 ( v_x597 | ( v_x597 << 2 ) )
#define v_x599 ( v_x595 | v_x598 )
#define v_x600 ( x & 3145728 )
#define v_x601 ( v_x600 << 20 )
#define v_x602 ( v_x601 | ( v_x601 << 2 ) )
#define v_x603 ( x & 12582912 )
#define v_x604 ( v_x603 << 22 )
#define v_x605 ( v_x604 | ( v_x604 << 2 ) )
#define v_x606 ( v_x602 | v_x605 )
#define v_x607 ( v_x599 | v_x606 )
#define v_x608 ( x & 50331648 )
#define v_x609 ( v_x608 << 24 )
#define v_x610 ( v_x609 | ( v_x609 << 2 ) )
#define v_x611 ( x & 201326592 )
#define v_x612 ( v_x611 << 26 )
#define v_x613 ( v_x612 | ( v_x612 << 2 ) )
#define v_x614 ( v_x610 | v_x613 )
#define v_x615 ( x & 805306368 )
#define v_x616 ( v_x615 << 28 )
#define v_x617 ( v_x616 | ( v_x616 << 2 ) )
#define v_x618 ( x & 3221225472 )
#define v_x619 ( v_x618 << 30 )
#define v_x620 ( v_x619 | ( v_x619 << 2 ) )
#define v_x621 ( v_x617 | v_x620 )
#define v_x622 ( v_x614 | v_x621 )
#define v_x623 ( v_x607 | v_x622 )
#define v_x624 ( v_x222 | v_x623 )
  case 125:
    return v_x624;
#define v_x625 ( x & 65536 )
#define v_x626 ( v_x625 << 16 )
#define v_x627 ( v_x626 | ( v_x626 << 1 ) )
#define v_x628 ( x & 131072 )
#define v_x629 ( v_x628 << 17 )
#define v_x630 ( v_x629 | ( v_x629 << 1 ) )
#define v_x631 ( v_x627 | v_x630 )
#define v_x632 ( x & 262144 )
#define v_x633 ( v_x632 << 18 )
#define v_x634 ( v_x633 | ( v_x633 << 1 ) )
#define v_x635 ( x & 524288 )
#define v_x636 ( v_x635 << 19 )
#define v_x637 ( v_x636 | ( v_x636 << 1 ) )
#define v_x638 ( v_x634 | v_x637 )
#define v_x639 ( v_x631 | v_x638 )
#define v_x640 ( x & 1048576 )
#define v_x641 ( v_x640 << 20 )
#define v_x642 ( v_x641 | ( v_x641 << 1 ) )
#define v_x643 ( x & 2097152 )
#define v_x644 ( v_x643 << 21 )
#define v_x645 ( v_x644 | ( v_x644 << 1 ) )
#define v_x646 ( v_x642 | v_x645 )
#define v_x647 ( x & 4194304 )
#define v_x648 ( v_x647 << 22 )
#define v_x649 ( v_x648 | ( v_x648 << 1 ) )
#define v_x650 ( x & 8388608 )
#define v_x651 ( v_x650 << 23 )
#define v_x652 ( v_x651 | ( v_x651 << 1 ) )
#define v_x653 ( v_x649 | v_x652 )
#define v_x654 ( v_x646 | v_x653 )
#define v_x655 ( v_x639 | v_x654 )
#define v_x656 ( x & 16777216 )
#define v_x657 ( v_x656 << 24 )
#define v_x658 ( v_x657 | ( v_x657 << 1 ) )
#define v_x659 ( x & 33554432 )
#define v_x660 ( v_x659 << 25 )
#define v_x661 ( v_x660 | ( v_x660 << 1 ) )
#define v_x662 ( v_x658 | v_x661 )
#define v_x663 ( x & 67108864 )
#define v_x664 ( v_x663 << 26 )
#define v_x665 ( v_x664 | ( v_x664 << 1 ) )
#define v_x666 ( x & 134217728 )
#define v_x667 ( v_x666 << 27 )
#define v_x668 ( v_x667 | ( v_x667 << 1 ) )
#define v_x669 ( v_x665 | v_x668 )
#define v_x670 ( v_x662 | v_x669 )
#define v_x671 ( x & 268435456 )
#define v_x672 ( v_x671 << 28 )
#define v_x673 ( v_x672 | ( v_x672 << 1 ) )
#define v_x674 ( x & 536870912 )
#define v_x675 ( v_x674 << 29 )
#define v_x676 ( v_x675 | ( v_x675 << 1 ) )
#define v_x677 ( v_x673 | v_x676 )
#define v_x678 ( x & 1073741824 )
#define v_x679 ( v_x678 << 30 )
#define v_x680 ( v_x679 | ( v_x679 << 1 ) )
#define v_x681 ( x & 2147483648 )
#define v_x682 ( v_x681 << 31 )
#define v_x683 ( v_x682 | ( v_x682 << 1 ) )
#define v_x684 ( v_x680 | v_x683 )
#define v_x685 ( v_x677 | v_x684 )
#define v_x686 ( v_x670 | v_x685 )
#define v_x687 ( v_x655 | v_x686 )
#define v_x688 ( v_x254 | v_x687 )
  case 126:
    return v_x688;
  case 127:
    return x;
  default:
    assert( false );
    return 0;
  }
}

} // namespace detail
} // namespace kitty

/*! \endcond */
