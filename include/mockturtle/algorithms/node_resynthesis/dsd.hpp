/* mockturtle: C++ logic network library
 * Copyright (C) 2018-2019  EPFL
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
  \file dsd.hpp
  \brief Use DSD as pre-process to resynthesis

  \author Mathias Soeken
*/

#pragma once

#include <vector>

#include <kitty/dynamic_truth_table.hpp>

#include "../dsd_decomposition.hpp"

namespace mockturtle
{

template<class Ntk, class ResynthesisFn>
class dsd_resynthesis
{
public:
  explicit dsd_resynthesis( ResynthesisFn& resyn_fn )
      : _resyn_fn( resyn_fn )
  {
  }

  template<typename LeavesIterator, typename Fn>
  void operator()( Ntk& ntk, kitty::dynamic_truth_table const& function, LeavesIterator begin, LeavesIterator end, Fn&& fn )
  {
    bool success{false};
    const auto on_prime = [&]( kitty::dynamic_truth_table const& remainder, std::vector<signal<Ntk>> const& leaves ) {
      signal<Ntk> f = ntk.get_constant( false );

      const auto on_signal = [&]( signal<Ntk> const& _f ) {
        if ( !success )
        {
          f = _f;
          success = true;
        }
        return true;
      };
      auto _leaves = leaves;
      _resyn_fn( ntk, remainder, _leaves.begin(), _leaves.end(), on_signal );
      return f;
    };

    const auto f = dsd_decomposition( ntk, function, std::vector<signal<Ntk>>( begin, end ), on_prime );
    if ( success )
    {
      fn( f );
    }
  }

private:
  ResynthesisFn& _resyn_fn;
};

//template<class Ntk, class Fn>
//dsd_resynthesis(Fn&&) -> dsd_resynthesis<Ntk, Fn>;

} /* namespace mockturtle */
