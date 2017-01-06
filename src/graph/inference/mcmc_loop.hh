// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2016 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef MCMC_LOOP_HH
#define MCMC_LOOP_HH

#include "config.h"

#include <iostream>
#include <queue>

#include <tuple>

#include "hash_map_wrap.hh"
#include "parallel_rng.hh"

#ifdef USING_OPENMP
#include <omp.h>
#endif
namespace graph_tool
{

template <class MCMCState, class RNG>
auto mcmc_sweep(MCMCState state, RNG& rng_)
{
    auto& g = state._g;

    vector<std::shared_ptr<RNG>> rngs;
    std::vector<std::pair<size_t, double>> best_move;

    if (state._parallel)
    {
        init_rngs(rngs, rng_);
        init_cache(state._E);
        best_move.resize(num_vertices(g));
    }

    auto& vlist = state._vlist;
    auto& beta = state._beta;


    double S = 0;
    size_t nmoves = 0;

    for (size_t iter = 0; iter < state._niter; ++iter)
    {
        if (!state._parallel)
        {
            std::shuffle(vlist.begin(), vlist.end(), rng_);
        }
        else
        {
            parallel_loop(vlist,
                          [&](size_t, auto v)
                          {
                              best_move[v] =
                                  std::make_pair(state.node_state(v),
                                                 numeric_limits<double>::max());
                          });
        }

        #pragma omp parallel firstprivate(state) if (state._parallel)
        parallel_loop_no_spawn
            (vlist,
             [&](size_t, auto v)
             {
                 auto& rng = get_rng(rngs, rng_);

                 if (!state._sequential)
                     v = uniform_sample(vlist, rng);

                 if (state.node_weight(v) == 0)
                     return;

                 auto r = state.node_state(v);
                 auto s = state.move_proposal(v, rng);

                 if (s == r)
                     return;

                 std::pair<double, double> dS = state.virtual_move_dS(v, s);

                 bool accept = false;
                 if (std::isinf(beta))
                 {
                     accept = dS.first < 0;
                 }
                 else
                 {
                     double a = -dS.first * beta + dS.second;
                     if (a > 0)
                     {
                         accept = true;
                     }
                     else
                     {
                         typedef std::uniform_real_distribution<> rdist_t;
                         double sample = rdist_t()(rng);
                         accept = sample < exp(a);
                     }
                 }

                 if (accept)
                 {
                     if (!state._parallel)
                     {
                         state.perform_move(v, s);
                         nmoves += state.node_weight(v);
                         S += dS.first;
                     }
                     else
                     {
                         best_move[v].first = s;
                         best_move[v].second = dS.first;
                     }
                     if (state._verbose)
                         cout << v << ": " << r << " -> " << s << " " << S << endl;
                 }
             });

        if (state._parallel)
        {
            for (auto v : vlist)
            {
                auto s = best_move[v].first;
                double dS = best_move[v].second;
                if (dS != numeric_limits<double>::max())
                {
                    dS = state.virtual_move_dS(v, s).first;

                    if (dS > 0 && std::isinf(beta))
                        continue;

                    state.perform_move(v, s);
                    nmoves++;
                    S += dS;
                }
            }
        }
    }
    return make_pair(S, nmoves);
}

} // graph_tool namespace

#endif //MCMC_LOOP_HH
