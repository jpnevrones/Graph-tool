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

#ifndef GRAPH_BLOCKMODEL_UTIL_HH
#define GRAPH_BLOCKMODEL_UTIL_HH

#include "config.h"

#include <tuple>

#include "hash_map_wrap.hh"

#include "../generation/sampler.hh"
#include "../generation/dynamic_sampler.hh"

#include "util.hh"
#include "int_part.hh"
#include "cache.hh"

#include <boost/multi_array.hpp>
#include <boost/math/special_functions/gamma.hpp>

#ifdef USING_OPENMP
#include <omp.h>
#endif

namespace graph_tool
{

// tuple utils
template <size_t i, class T>
void add_to_tuple_imp(T&)
{
}

template <size_t i, class T, class Ti, class... Ts>
void add_to_tuple_imp(T& tuple, Ti&& v, Ts&&... vals)
{
    get<i>(tuple) += v;
    add_to_tuple_imp<i+1>(tuple, vals...);
}

template <class T, class... Ts>
void add_to_tuple(T& tuple, Ts&&... vals)
{
    add_to_tuple_imp<0>(tuple, vals...);
}

namespace detail {
   template <class F, class Tuple, std::size_t... I>
   constexpr decltype(auto) tuple_apply_impl(F&& f, Tuple&& t,
                                             std::index_sequence<I...>)
   {
       return f(std::get<I>(std::forward<Tuple>(t))...);
   }
} // namespace detail

template <class F, class Tuple>
constexpr decltype(auto) tuple_apply(F&& f, Tuple&& t)
{
    return detail::tuple_apply_impl
        (std::forward<F>(f), std::forward<Tuple>(t),
         std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>{}>{});
}

// ====================
// Entropy calculation
// ====================

enum deg_dl_kind
{
    ENT,
    UNIFORM,
    DIST
};

struct entropy_args_t
{
    bool dense;
    bool multigraph;
    bool exact;
    bool adjacency;
    bool partition_dl;
    bool degree_dl;
    deg_dl_kind  degree_dl_kind;
    bool edges_dl;
};

// Sparse entropy terms
// ====================

// exact microcanonical deg-corr entropy
template <class Graph>
inline double eterm_exact(size_t r, size_t s, size_t mrs, const Graph&)
{
    double val = lgamma_fast(mrs + 1);

    if (is_directed::apply<Graph>::type::value || r != s)
    {
        return -val;
    }
    else
    {
#ifndef __clang__
        constexpr double log_2 = log(2);
#else
        const double log_2 = log(2);
#endif
        return -val - mrs * log_2;
    }
}

template <class Graph>
inline double vterm_exact(size_t mrp, size_t mrm, size_t wr, bool deg_corr,
                          const Graph&)
{
    if (deg_corr)
    {
        if (is_directed::apply<Graph>::type::value)
            return lgamma_fast(mrp + 1) + lgamma_fast(mrm + 1);
        else
            return lgamma_fast(mrp + 1);
    }
    else
    {
        if (is_directed::apply<Graph>::type::value)
            return (mrp + mrm) * safelog(wr);
        else
            return mrp * safelog(wr);
    }
}

// "edge" term of the entropy
template <class Graph>
inline double eterm(size_t r, size_t s, size_t mrs, const Graph&)
{
    if (!is_directed::apply<Graph>::type::value && r == s)
        mrs *= 2;

    double val = xlogx(mrs);

    if (is_directed::apply<Graph>::type::value || r != s)
        return -val;
    else
        return -val / 2;
}

// "vertex" term of the entropy
template <class Graph>
inline double vterm(size_t mrp, size_t mrm, size_t wr, bool deg_corr,
                    Graph&)
{
    double one = 0.5;

    if (is_directed::apply<Graph>::type::value)
        one = 1;

    if (deg_corr)
        return one * (xlogx(mrm) + xlogx(mrp));
    else
        return one * (mrm * safelog(wr) + mrp * safelog(wr));
}



// Dense entropy
// =============

// "edge" term of the entropy
template <class Graph>
inline double eterm_dense(size_t r, size_t s, int ers, double wr_r,
                          double wr_s, bool multigraph, const Graph&)
{
    // we should not use integers here, since they may overflow
    double nrns;

    if (ers == 0)
        return 0.;

    assert(wr_r + wr_s > 0);

    if (r != s || is_directed::apply<Graph>::type::value)
    {
        nrns = wr_r * wr_s;
    }
    else
    {
        if (multigraph)
            nrns = (wr_r * (wr_r + 1)) / 2;
        else
            nrns = (wr_r * (wr_r - 1)) / 2;
    }

    double S;
    if (multigraph)
        S = lbinom(nrns + ers - 1, ers); // do not use lbinom_fast!
    else
        S = lbinom(nrns, ers);
    return S;
}

// Weighted entropy terms

template <class DT>
double positive_w_log_P(DT N, double x, double alpha, double beta)
{
    if (N == 0)
        return 0.;
    return boost::math::lgamma(N + alpha) - boost::math::lgamma(alpha)
        + alpha * log(beta) - (alpha + N) * log(beta + x);
}

template <class DT>
double signed_w_log_P(DT N, double x, double v, double m0, double k0, double v0,
                      double nu0)
{
    if (N == 0)
        return 0.;
    auto k_n = k0 + N;
    auto nu_n = nu0 + N;
    auto v_n = (v0 * nu0 + v + ((N * k0)/(k0 + N)) * pow(m0 - x/N, 2.)) / nu_n;
    return boost::math::lgamma(nu_n/2.) - boost::math::lgamma(nu0/2.)
        + (log(k0) - log(k_n))/2. + (nu0 / 2.) * log(nu0 * v0)
        - (nu_n / 2.) * log(nu_n * v_n) - (N/2.) * log(M_PI);
}


// ===============
// Partition stats
// ===============

constexpr size_t null_group = std::numeric_limits<size_t>::max();

typedef vprop_map_t<std::vector<std::tuple<size_t, size_t, size_t>>>::type
    degs_map_t;

struct simple_degs_t {};

template <class Graph, class Vprop, class Eprop>
auto get_degs(size_t v, Vprop& vweight, Eprop& eweight, const simple_degs_t&,
              Graph& g)
{
    std::array<std::tuple<size_t, size_t, size_t>, 1>
        degs {{make_tuple(size_t(in_degreeS()(v, g, eweight)),
                          size_t(out_degreeS()(v, g, eweight)),
                          size_t(vweight[v]))}};
    return degs;
}

template <class Graph, class Vprop, class Eprop>
auto& get_degs(size_t v, Vprop&, Eprop&, typename degs_map_t::unchecked_t& degs,
               Graph&)
{
    return degs[v];
}

class partition_stats_t
{
public:

    typedef gt_hash_map<pair<size_t,size_t>, int> map_t;

    template <class Graph, class Vprop, class VWprop, class Eprop, class Degs,
              class Mprop, class Vlist>
    partition_stats_t(Graph& g, Vprop& b, Vlist& vlist, size_t E, size_t B,
                      VWprop& vweight, Eprop& eweight, Degs& degs,
                      const Mprop& ignore_degree, std::vector<size_t>& bmap,
                      bool allow_empty)
        : _bmap(bmap), _N(0), _E(E), _total_B(B), _allow_empty(allow_empty)
    {
        for (auto v : vlist)
        {
            if (vweight[v] == 0)
                continue;

            auto r = get_r(b[v]);

            auto&& ks = get_degs(v, vweight, eweight, degs, g);
            if (v >= _ignore_degree.size())
                _ignore_degree.resize(v + 1, 0);
            _ignore_degree[v] = ignore_degree[v];
            for (auto& k : ks)
            {
                size_t kin = get<0>(k);
                size_t kout = get<1>(k);
                size_t n = get<2>(k);
                if (_ignore_degree[v] == 2)
                    kout = 0;
                if (_ignore_degree[v] != 1)
                {
                    _hist[r][make_pair(kin, kout)] += n;
                    _em[r] += kin * n;
                    _ep[r] += kout * n;
                }
                _total[r] += n;
                _N += n;
            }
        }

        _actual_B = 0;
        for (auto n : _total)
        {
            if (n > 0)
                _actual_B++;
        }
    }

    size_t get_r(size_t r)
    {
        constexpr size_t null =
            std::numeric_limits<size_t>::max();
        if (r >= _bmap.size())
            _bmap.resize(r + 1, null);
        size_t nr = _bmap[r];
        if (nr == null)
            nr = _bmap[r] = _hist.size();
        if (nr >= _hist.size())
        {
            _hist.resize(nr + 1);
            _total.resize(nr + 1);
            _ep.resize(nr + 1);
            _em.resize(nr + 1);
        }
        return nr;
    }

    double get_partition_dl()
    {
        double S = 0;
        if (_allow_empty)
            S += lbinom(_total_B + _N - 1, _N);
        else
            S += lbinom(_N - 1, _actual_B - 1);
        S += lgamma_fast(_N + 1);
        for (auto nr : _total)
            S -= lgamma_fast(nr + 1);
        return S;
    }

    double get_deg_dl_ent()
    {
        double S = 0;
        for (size_t r = 0; r < _ep.size(); ++r)
        {
            size_t total = 0;
            for (auto& k_c : _hist[r])
            {
                S -= xlogx(k_c.second);
                total += k_c.second;
            }
            S += xlogx(total);
        }
        return S;
    }

    double get_deg_dl_uniform()
    {
        double S = 0;
        for (size_t r = 0; r < _ep.size(); ++r)
        {
            S += lbinom(_total[r] + _ep[r] - 1, _ep[r]);
            S += lbinom(_total[r] + _em[r] - 1, _em[r]);
        }
        return S;
    }

    double get_deg_dl_dist()
    {
        double S = 0;
        for (size_t r = 0; r < _ep.size(); ++r)
        {
            S += log_q(_ep[r], _total[r]);
            S += log_q(_em[r], _total[r]);

            size_t total = 0;
            for (auto& k_c : _hist[r])
            {
                S -= lgamma_fast(k_c.second + 1);
                total += k_c.second;
            }
            S += lgamma_fast(total + 1);
        }
        return S;
    }

    double get_deg_dl(int kind)
    {
        switch (kind)
        {
        case deg_dl_kind::ENT:
            return get_deg_dl_ent();
        case deg_dl_kind::UNIFORM:
            return get_deg_dl_uniform();
        case deg_dl_kind::DIST:
            return get_deg_dl_dist();
        default:
            return numeric_limits<double>::quiet_NaN();
        }
    }

    template <class VProp>
    double get_delta_partition_dl(size_t v, size_t r, size_t nr, VProp& vweight)
    {
        if (r == nr)
            return 0;

        if (r != null_group)
            r = get_r(r);

        if (nr != null_group)
            nr = get_r(nr);

        int n = vweight[v];
        if (n == 0)
        {
            if (r == null_group)
                n = 1;
            else
                return 0;
        }

        double S_b = 0, S_a = 0;

        if (r != null_group)
        {
            S_b += -lgamma_fast(_total[r] + 1);
            S_a += -lgamma_fast(_total[r] - n + 1);
        }

        if (nr != null_group)
        {
            S_b += -lgamma_fast(_total[nr] + 1);
            S_a += -lgamma_fast(_total[nr] + n + 1);
        }

        int dN = 0;
        if (r == null_group)
            dN += n;
        if (nr == null_group)
            dN -= n;

        S_b += lgamma_fast(_N + 1);
        S_a += lgamma_fast(_N + dN + 1);

        int dB = 0;
        if (r != null_group && _total[r] == n)
            dB--;
        if (nr != null_group && _total[nr] == 0)
            dB++;

        if ((dN != 0 || dB != 0) && !_allow_empty)
        {
            //S_b += lbinom_fast(std::min(_total_B, _N), _actual_B);
            S_b += lbinom_fast(_N - 1, _actual_B - 1);
            //S_a += lbinom_fast(std::min(_total_B, _N + dN), _actual_B + dB);
            S_a += lbinom_fast(_N - 1 + dN, _actual_B + dB - 1);
        }

        return S_a - S_b;
    }

    template <class Graph, class VProp>
    double get_delta_edges_dl(size_t v, size_t r, size_t nr, VProp& vweight,
                              Graph&)
    {
        if (r == nr || _allow_empty)
            return 0;

        if (r != null_group)
            r = get_r(r);
        if (nr != null_group)
            nr = get_r(nr);

        double S_b = 0, S_a = 0;

        int n = vweight[v];

        if (n == 0)
        {
            if (r == null_group)
                n = 1;
            else
                return 0;
        }

        int dB = 0;
        if (r != null_group && _total[r] == n)
            dB--;
        if (nr != null_group && _total[nr] == 0)
            dB++;

        if (dB != 0)
        {
            auto get_x = [](size_t B)
                {
                    if (is_directed::apply<Graph>::type::value)
                        return B * B;
                    else
                        return (B * (B + 1)) / 2;
                };

            S_b += lbinom(get_x(_actual_B) + _E - 1, _E);
            S_a += lbinom(get_x(_actual_B + dB) + _E - 1, _E);
        }

        return S_a - S_b;
    }

    template <class Graph, class VProp, class EProp, class Degs>
    double get_delta_deg_dl(size_t v, size_t r, size_t nr, VProp& vweight,
                            EProp& eweight, Degs& degs, Graph& g, int kind)
    {
        if (r == nr || _ignore_degree[v] == 1 || vweight[v] == 0)
            return 0;
        if (r != null_group)
            r = get_r(r);
        if (nr != null_group)
            nr = get_r(nr);
        auto&& ks = get_degs(v, vweight, eweight, degs, g);
        auto* _ks = &ks;
        typename std::remove_reference<decltype(ks)>::type nks;
        if (_ignore_degree[v] == 2)
        {
            nks = ks;
            for (auto& k : nks)
                get<1>(k) = 0;
            _ks = &nks;
        }
        double dS = 0;
        switch (kind)
        {
        case deg_dl_kind::ENT:
            if (r != null_group)
                dS += get_delta_deg_dl_ent_change(r,  *_ks, -1);
            if (nr != null_group)
                dS += get_delta_deg_dl_ent_change(nr, *_ks, +1);
            break;
        case deg_dl_kind::UNIFORM:
            if (r != null_group)
                dS += get_delta_deg_dl_uniform_change(v, r,  *_ks, -1);
            if (nr != null_group)
                dS += get_delta_deg_dl_uniform_change(v, nr, *_ks, +1);
            break;
        case deg_dl_kind::DIST:
            if (r != null_group)
                dS += get_delta_deg_dl_dist_change(v, r,  *_ks, -1);
            if (nr != null_group)
                dS += get_delta_deg_dl_dist_change(v, nr, *_ks, +1);
            break;
        default:
            dS = numeric_limits<double>::quiet_NaN();
        }
        return dS;
    }

    template <class Ks>
    double get_delta_deg_dl_ent_change(size_t r, Ks& ks, int diff)
    {
        int nr = _total[r];
        auto get_Sk = [&](size_t s, pair<size_t, size_t>& deg, int delta)
            {
                int nd = 0;
                auto iter = _hist[s].find(deg);
                if (iter != _hist[s].end())
                    nd = iter->second;
                assert(nd + delta >= 0);
                return -xlogx(nd + delta);
            };

        double S_b = 0, S_a = 0;
        int dn = 0;
        for (auto& k : ks)
        {
            size_t kin = get<0>(k);
            size_t kout = get<1>(k);
            int nk = get<2>(k);
            dn += diff * nk;
            auto deg = make_pair(kin, kout);
            S_b += get_Sk(r, deg,         0);
            S_a += get_Sk(r, deg, diff * nk);
        }
        S_b += xlogx(nr);
        S_a += xlogx(nr + dn);

        return S_a - S_b;
    }

    template <class Ks>
    double get_delta_deg_dl_uniform_change(size_t v, size_t r, Ks& ks, int diff)
    {
        auto get_Se = [&](int dn, int dkin, int dkout)
            {
                double S = 0;
                S += lbinom(_total[r] + dn + _ep[r] - 1 + dkout, _ep[r] + dkout);
                S += lbinom(_total[r] + dn + _em[r] - 1 + dkin,  _em[r] + dkin);
                return S;
            };

        double S_b = 0, S_a = 0;
        int tkin = 0, tkout = 0, n = 0;
        for (auto& k : ks)
        {
            size_t kin = get<0>(k);
            size_t kout = get<1>(k);
            int nk = get<2>(k);
            tkin += kin * nk;
            if (_ignore_degree[v] != 2)
                tkout += kout * nk;
            n += nk;
        }

        S_b += get_Se(       0,           0,            0);
        S_a += get_Se(diff * n, diff * tkin, diff * tkout);
        return S_a - S_b;
    }


    template <class Ks>
    double get_delta_deg_dl_dist_change(size_t v, size_t r, Ks& ks, int diff)
    {

        auto get_Se = [&](int delta, int kin, int kout)
            {
                double S = 0;
                assert(_total[r] + delta >= 0);
                assert(_em[r] + kin >= 0);
                assert(_ep[r] + kout >= 0);
                S += log_q(_em[r] + kin, _total[r] + delta);
                S += log_q(_ep[r] + kout, _total[r] + delta);
                return S;
            };

        auto get_Sr = [&](int delta)
            {
                assert(_total[r] + delta + 1 >= 0);
                return lgamma_fast(_total[r] + delta + 1);
            };

        auto get_Sk = [&](pair<size_t, size_t>& deg, int delta)
            {
                int nd = 0;
                auto iter = _hist[r].find(deg);
                if (iter != _hist[r].end())
                    nd = iter->second;
                assert(nd + delta >= 0);
                return -lgamma_fast(nd + delta + 1);
            };

        double S_b = 0, S_a = 0;
        int tkin = 0, tkout = 0, n = 0;
        for (auto& k : ks)
        {
            size_t kin = get<0>(k);
            size_t kout = get<1>(k);
            int nk = get<2>(k);
            tkin += kin * nk;
            if (_ignore_degree[v] != 2)
                tkout += kout * nk;
            n += nk;

            auto deg = make_pair(kin, kout);
            S_b += get_Sk(deg,         0);
            S_a += get_Sk(deg, diff * nk);
        }

        S_b += get_Se(       0,           0,            0);
        S_a += get_Se(diff * n, diff * tkin, diff * tkout);

        S_b += get_Sr(       0);
        S_a += get_Sr(diff * n);

        return S_a - S_b;
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void change_vertex(size_t v, size_t r, bool deg_corr, Graph& g,
                       VWeight& vweight, EWeight& eweight, Degs& degs,
                       int diff)
    {
        auto&& ks = get_degs(v, vweight, eweight, degs, g);
        for (auto& k : ks)
        {
            auto kin = get<0>(k);
            auto kout = get<1>(k);
            int n = get<2>(k);
            change_k(v, r, deg_corr, n, kin, kout, diff);
        }
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void remove_vertex(size_t v, size_t r, bool deg_corr, Graph& g,
                       VWeight& vweight, EWeight& eweight, Degs& degs)
    {
        if (r == null_group || vweight[v] == 0)
            return;
        r = get_r(r);
        change_vertex(v, r, deg_corr, g, vweight, eweight, degs, -1);
    }

    template <class Graph, class VWeight, class EWeight, class Degs>
    void add_vertex(size_t v, size_t nr, bool deg_corr, Graph& g,
                    VWeight& vweight, EWeight& eweight, Degs& degs)
    {
        if (nr == null_group || vweight[v] == 0)
            return;
        nr = get_r(nr);
        change_vertex(v, nr, deg_corr, g, vweight, eweight, degs, 1);
    }

    void change_k(size_t v, size_t r, bool deg_corr, int vweight,
                  int kin, int kout, int diff)
    {
        if (_total[r] == 0 && diff * vweight > 0)
            _actual_B++;

        if (_total[r] == vweight && diff * vweight < 0)
            _actual_B--;

        _total[r] += diff * vweight;
        _N += diff * vweight;

        assert(_total[r] >= 0);

        if (deg_corr && _ignore_degree[v] != 1)
        {
            if (_ignore_degree[v] == 2)
                kout = 0;
            auto deg = make_pair(kin, kout);
            _hist[r][deg] += diff * vweight;
            _em[r] += diff * deg.first * vweight;
            _ep[r] += diff * deg.second * vweight;
        }
    }

private:
    vector<size_t>& _bmap;
    size_t _N;
    size_t _E;
    size_t _actual_B;
    size_t _total_B;
    bool _allow_empty;
    vector<map_t> _hist;
    vector<int> _total;
    vector<int> _ep;
    vector<int> _em;
    vector<uint8_t> _ignore_degree;
};

// ===============================
// Block moves
// ===============================

// this structure speeds up the access to the edges between given blocks, since
// we're using an adjacency list to store the block structure (it is simply an
// adjacency matrix)

template <class Graph, class BGraph>
class EMat
{
public:
    template <class Vprop, class RNG>
    EMat(Graph& g, Vprop b, BGraph& bg, RNG&)
        : _bedge(get(edge_index_t(), g), 0)
    {
        sync(g, b, bg);
    }

    template <class Vprop>
    void sync(Graph& g, Vprop b, BGraph& bg)
    {
        size_t B = num_vertices(bg);
        _mat.resize(boost::extents[B][B]);
        std::fill(_mat.data(), _mat.data() + _mat.num_elements(), _null_edge);

        for (auto e : edges_range(bg))
        {
            assert(get_me(source(e, bg),target(e, bg)) == _null_edge);
            _mat[source(e, bg)][target(e, bg)] = e;
            if (!is_directed::apply<BGraph>::type::value)
                _mat[target(e, bg)][source(e, bg)] = e;
        }

        auto bedge_c = _bedge.get_checked();
        for (auto e : edges_range(g))
        {
            auto r = b[source(e, g)];
            auto s = b[target(e, g)];
            bedge_c[e] = _mat[r][s];
            assert(bedge_c[e] != _null_edge);
        }
    }

    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<BGraph>::edge_descriptor edge_t;

    const auto& get_me(vertex_t r, vertex_t s) const
    {
        return _mat[r][s];
    }

    void put_me(vertex_t r, vertex_t s, const edge_t& e)
    {
        _mat[r][s] = e;
        if (!is_directed::apply<BGraph>::type::value && r != s)
            _mat[s][r] = e;
    }

    void remove_me(vertex_t r, vertex_t s, const edge_t& me, BGraph& bg,
                   bool delete_edge = true)
    {
        if (delete_edge)
        {
            _mat[r][s] = _null_edge;
            if (!is_directed::apply<Graph>::type::value)
                _mat[s][r] = _null_edge;
            remove_edge(me, bg);
        }
    }

    template <class Edge>
    auto& get_bedge(const Edge& e) { return _bedge[e]; }
    auto& get_bedge_map() { return _bedge; }
    const auto& get_bedge_map() const { return _bedge; }
    const auto& get_null_edge() const { return _null_edge; }

private:
    multi_array<edge_t, 2> _mat;
    typedef typename property_map_type::apply
        <edge_t,
         typename property_map<Graph,
                               edge_index_t>::type>::type bedge_t;
    typename bedge_t::unchecked_t _bedge;
    static const edge_t _null_edge;
};

template <class Graph, class BGraph>
const typename EMat<Graph, BGraph>::edge_t EMat<Graph, BGraph>::_null_edge;


template <class Key>
class perfect_hash_t
{
public:
    template <class RNG>
    perfect_hash_t(size_t N, RNG& rng)
        : _index(std::make_shared<std::vector<size_t>>())
    {
        auto& index = *_index;
        index.reserve(N);
        for (size_t i = 0; i < N; ++i)
            index.push_back(i);
        std::shuffle(index.begin(), index.end(), rng);
    }
    perfect_hash_t() {}
    size_t operator()(const Key& k) const {return (*_index)[k];}
private:
    std::shared_ptr<std::vector<size_t>> _index;
};

// this structure speeds up the access to the edges between given blocks, since
// we're using an adjacency list to store the block structure (this is like
// EMat above, but takes less space and is slower)

template <class Graph, class BGraph>
class EHash
{
public:

    template <class Vprop, class RNG>
    EHash(Graph& g, Vprop b, BGraph& bg, RNG& rng)
        : _hash_function(num_vertices(bg), rng),
          _hash(num_vertices(bg), ehash_t(0, _hash_function)),
          _bedge(get(edge_index_t(), g), 0)
    {
        sync(g, b, bg);
    }

    template <class Vprop>
    void sync(Graph& g, Vprop b, BGraph& bg)
    {
        _hash.clear();
        _hash.resize(num_vertices(bg), ehash_t(0, _hash_function));

        for (auto e : edges_range(bg))
        {
            assert(get_me(source(e, bg), target(e, bg)) == _null_edge);
            put_me(source(e, bg), target(e, bg), e);
        }

        auto bedge_c = _bedge.get_checked();
        for (auto e : edges_range(g))
        {
            auto r = b[source(e, g)];
            auto s = b[target(e, g)];
            bedge_c[e] = get_me(r, s);
            assert(bedge_c[e] != _null_edge);
        }
    }

    typedef typename graph_traits<BGraph>::vertex_descriptor vertex_t;
    typedef typename graph_traits<BGraph>::edge_descriptor edge_t;

    const auto& get_me(vertex_t r, vertex_t s) const
    {
        auto& map = _hash[r];
        const auto& iter = map.find(s);
        if (iter == map.end())
            return _null_edge;
        return iter->second;
    }

    void put_me(vertex_t r, vertex_t s, const edge_t& e)
    {
        assert(r < _hash.size());
        _hash[r][s] = e;
        if (!is_directed::apply<Graph>::type::value)
            _hash[s][r] = e;
    }

    void remove_me(vertex_t r, vertex_t s, const edge_t& me, BGraph& bg,
                   bool delete_edge = true)
    {
        if (delete_edge)
        {
            assert(r < _hash.size());
            _hash[r].erase(s);
            if (!is_directed::apply<Graph>::type::value)
                _hash[s].erase(r);
            remove_edge(me, bg);
        }
    }

    template <class Edge>
    auto& get_bedge(const Edge& e) { return _bedge[e]; }
    auto& get_bedge_map() { return _bedge; }
    const auto& get_bedge_map() const { return _bedge; }
    const auto& get_null_edge() const { return _null_edge; }

private:
    perfect_hash_t<vertex_t> _hash_function;
    typedef gt_hash_map<vertex_t, edge_t, perfect_hash_t<vertex_t>> ehash_t;
    std::vector<ehash_t> _hash;
    typedef typename property_map_type::apply
        <edge_t,
         typename property_map<Graph,
                               edge_index_t>::type>::type bedge_t;
    typename bedge_t::unchecked_t _bedge;
    static const edge_t _null_edge;
};

template <class Graph, class BGraph>
const typename EHash<Graph, BGraph>::edge_t EHash<Graph, BGraph>::_null_edge;

template <class Vertex, class Eprop, class Emat, class BEdge>
inline auto get_beprop(Vertex r, Vertex s, const Eprop& eprop, const Emat& emat,
                       BEdge& me)
{
    typedef typename property_traits<Eprop>::value_type val_t;
    me = emat.get_me(r, s);
    if (me != emat.get_null_edge())
        return eprop[me];
    return val_t();
}

template <class Vertex, class Eprop, class Emat>
inline auto get_beprop(Vertex r, Vertex s, const Eprop& eprop, const Emat& emat)
{
    typedef typename property_traits<Eprop>::key_type bedge_t;
    bedge_t me;
    return get_beprop(r, s, eprop, emat, me);
}


// Manage a set of block pairs and corresponding edge counts that will be
// updated

template <class Graph, class BGraph, class... EVals>
class EntrySet
{
public:
    typedef typename graph_traits<BGraph>::edge_descriptor bedge_t;

    EntrySet(size_t B)
    {
        _r_field_t.resize(B, _null);
        _nr_field_t.resize(B, _null);

        if (is_directed::apply<Graph>::type::value)
        {
            _r_field_s.resize(B, _null);
            _nr_field_s.resize(B, _null);
        }
    }

    void set_move(size_t r, size_t nr)
    {
        _rnr = make_pair(r, nr);
    }

    template <class... DVals>
    void insert_delta(size_t t, size_t s, const bedge_t& me, DVals... delta)
    {
        insert_delta_imp(t, s, me, typename is_directed::apply<Graph>::type(),
                         delta...);
    }

    template <class... DVals>
    void insert_delta_imp(size_t t, size_t s, const bedge_t& me, std::true_type,
                          DVals... delta)
    {
        bool src = false;
        if (t != _rnr.first && t != _rnr.second)
        {
            std::swap(t, s);
            src = true;
        }

        assert(t == _rnr.first || t == _rnr.second);

        auto& r_field = (src) ? _r_field_s : _r_field_t;
        auto& nr_field = (src) ? _nr_field_s : _nr_field_t;

        auto& field = (_rnr.first == t) ? r_field : nr_field;
        if (field[s] == _null)
        {
            field[s] = _entries.size();
            if (src)
                _entries.emplace_back(s, t);
            else
                _entries.emplace_back(t, s);
            _delta.emplace_back();
            _mes.push_back(me);
        }
        add_to_tuple(_delta[field[s]], delta...);
    }

    template <class... DVals>
    void insert_delta_imp(size_t t, size_t s, const bedge_t& me, std::false_type,
                          DVals... delta)
    {
        if (t > s)
            std::swap(t, s);
        if (t != _rnr.first && t != _rnr.second)
            std::swap(t, s);

        assert(t == _rnr.first || t == _rnr.second);

        auto& r_field = _r_field_t;
        auto& nr_field = _nr_field_t;

        auto& field = (_rnr.first == t) ? r_field : nr_field;
        if (field[s] == _null)
        {
            field[s] = _entries.size();
            _entries.emplace_back(t, s);
            _delta.emplace_back();
            _mes.push_back(me);
        }
        add_to_tuple(_delta[field[s]], delta...);
    }

    const auto& get_delta(size_t t, size_t s)
    {
        if (is_directed::apply<Graph>::type::value)
        {
            if (t == _rnr.first || t == _rnr.second)
                return get_delta_target(t, s);
            if (s == _rnr.first || s == _rnr.second)
                return get_delta_source(t, s);
            return _null_delta;
        }
        else
        {
            if (t > s)
                std::swap(t, s);
            if (t != _rnr.first && t != _rnr.second)
                std::swap(t, s);
            if (t == _rnr.first || t == _rnr.second)
                return get_delta_target(t, s);
            return _null_delta;
        }
    }

    const auto& get_delta_target(size_t r, size_t s)
    {
        vector<size_t>& field = (_rnr.first == r) ? _r_field_t : _nr_field_t;
        if (field[s] == _null)
            return _null_delta;
        else
            return _delta[field[s]];
    }

    const auto& get_delta_source(size_t s, size_t r)
    {
        vector<size_t>& field = (_rnr.first == r) ? _r_field_s : _nr_field_s;
        if (field[s] == _null)
            return _null_delta;
        else
            return _delta[field[s]];
    }

    void clear()
    {
        for (const auto& rs : _entries)
        {
            size_t r = rs.first;
            size_t s = rs.second;
            _r_field_t[r] = _nr_field_t[r] = _null;
            _r_field_t[s] = _nr_field_t[s] = _null;
            if (is_directed::apply<Graph>::type::value)
            {
                _r_field_s[r] = _nr_field_s[r] = _null;
                _r_field_s[s] = _nr_field_s[s] = _null;
            }
        }
        _entries.clear();
        _delta.clear();
        _mes.clear();
    }

    const vector<pair<size_t, size_t> >& get_entries() { return _entries; }
    const vector<std::tuple<EVals...>>& get_delta() { return _delta; }
    vector<bedge_t>& get_mes() { return _mes; }

    const auto& get_null_edge() const { return _null_edge; }

private:
    static constexpr size_t _null = numeric_limits<size_t>::max();
    static const std::tuple<EVals...> _null_delta;
    static const bedge_t _null_edge;

    pair<size_t, size_t> _rnr;
    vector<size_t> _r_field_t;
    vector<size_t> _nr_field_t;
    vector<size_t> _r_field_s;
    vector<size_t> _nr_field_s;
    vector<pair<size_t, size_t>> _entries;
    vector<std::tuple<EVals...>> _delta;
    vector<bedge_t> _mes;
};

template <class Graph, class BGraph, class... EVals>
constexpr size_t EntrySet<Graph, BGraph, EVals...>::_null;

template <class Graph, class BGraph, class... EVals>
const std::tuple<EVals...> EntrySet<Graph, BGraph, EVals...>::_null_delta;

template <class Graph, class BGraph, class... EVals>
const typename graph_traits<BGraph>::edge_descriptor
EntrySet<Graph, BGraph, EVals...>::_null_edge;

struct is_loop_nop
{
    bool operator()(size_t) const { return false; }
};

template <bool Add, class Graph, class BGraph, class Vertex, class Vprop,
          class Eprop, class EBedge, class MEntries, class IL,
          class... Eprops>
void modify_entries(Vertex v, Vertex r, Vprop& b, EBedge& bedge, Graph& g,
                    BGraph&, MEntries& m_entries, const IL& is_loop,
                    Eprop& eweights, Eprops&... eprops)
{
    typedef typename graph_traits<Graph>::vertex_descriptor vertex_t;
    std::tuple<int, typename property_traits<Eprops>::value_type...>
        self_weight;
    for (auto e : out_edges_range(v, g))
    {
        vertex_t u = target(e, g);
        vertex_t s = b[u];
        int ew = eweights[e];
        //assert(ew > 0);

        if (Add)
        {
            if (u == v)
                s = r;
            m_entries.insert_delta(r, s, m_entries.get_null_edge(),
                                   ew, eprops[e]...);
        }
        else
        {
            const auto& me = bedge[e];
            m_entries.insert_delta(r, s, me, -ew, -eprops[e]...);
        }

        if ((u == v || is_loop(v)) && !is_directed::apply<Graph>::type::value)
            add_to_tuple(self_weight, ew, eprops[e]...);
    }

    if (get<0>(self_weight) > 0 && get<0>(self_weight) % 2 == 0 &&
        !is_directed::apply<Graph>::type::value)
        tuple_apply([&](auto&&... vals)
                    {
                        if (Add)
                            m_entries.insert_delta(r, r,
                                                   m_entries.get_null_edge(),
                                                   (-vals / 2)...);
                        else
                            m_entries.insert_delta(r, r,
                                                   m_entries.get_null_edge(),
                                                   (vals / 2)...);
                    }, self_weight);

    for (auto e : in_edges_range(v, g))
    {
        vertex_t u = source(e, g);
        if (u == v)
            continue;
        vertex_t s = b[u];
        int ew = eweights[e];

        if (Add)
        {
            m_entries.insert_delta(s, r, m_entries.get_null_edge(),
                                   ew, eprops[e]...);
        }
        else
        {
            const auto& me = bedge[e];
            m_entries.insert_delta(s, r, me, -ew, -eprops[e]...);
        }
    }
}

// obtain the necessary entries in the e_rs matrix which need to be modified
// after the move
template <class Graph, class BGraph, class Vertex, class Vprop, class Eprop,
          class EBedge, class MEntries, class IL, class... Eprops>
void move_entries(Vertex v, size_t r, size_t nr, Vprop& b, EBedge& bedge,
                  Graph& g, BGraph& bg, MEntries& m_entries, const IL& is_loop,
                  Eprop& eweights, Eprops&... eprops)
{
    m_entries.set_move(r, nr);

    if (r != null_group)
        modify_entries<false>(v, r, b, bedge, g, bg, m_entries, is_loop,
                              eweights, eprops...);
    if (nr != null_group)
        modify_entries<true>(v, nr, b, bedge, g, bg, m_entries, is_loop,
                             eweights, eprops...);
}


// operation on a set of entries
template <class MEntries,  class OP>
void entries_op(MEntries& m_entries, OP&& op)
{
    const auto& entries = m_entries.get_entries();
    const auto& delta = m_entries.get_delta();
    auto& mes = m_entries.get_mes();

    for (size_t i = 0; i < entries.size(); ++i)
    {
        auto& entry = entries[i];
        auto er = entry.first;
        auto es = entry.second;
        auto& me = mes[i];
        op(er, es, me, delta[i]);
    }
}

// obtain the entropy difference given a set of entries in the e_rs matrix
template <bool exact, class MEntries, class Eprop, class EMat, class BGraph>
double entries_dS(MEntries& m_entries, Eprop& mrs, EMat& emat, BGraph& bg)
{
    double dS = 0;
    entries_op(m_entries,
               [&](auto r, auto s, auto& me, auto& delta)
               {
                   size_t ers;
                   if (me == m_entries.get_null_edge())
                       ers = get_beprop(r, s, mrs, emat, me); // slower
                   else
                       ers = mrs[me];
                   int d = get<0>(delta);
                   assert(int(ers) + d >= 0);
                   if (exact)
                       dS += eterm_exact(r, s, ers + d, bg) - eterm_exact(r, s, ers, bg);
                   else
                       dS += eterm(r, s, ers + d, bg) - eterm(r, s, ers, bg);
               });
    return dS;
}


// ====================================
// Construct and manage half-edge lists
// ====================================

//the following guarantees a stable (source, target) ordering even for
//undirected graphs
template <class Edge, class Graph>
inline typename graph_traits<Graph>::vertex_descriptor
get_source(const Edge& e, const Graph &g)
{
    if (is_directed::apply<Graph>::type::value)
        return source(e, g);
    return min(source(e, g), target(e, g));
}

template <class Edge, class Graph>
inline typename graph_traits<Graph>::vertex_descriptor
get_target(const Edge& e, const Graph &g)
{
    if (is_directed::apply<Graph>::type::value)
        return target(e, g);
    return max(source(e, g), target(e, g));
}


template <class Graph, class Weighted>
class EGroups
{
public:
    template <class Vprop, class Eprop, class BGraph>
    void init(Vprop b, Eprop eweight, Graph& g, BGraph& bg)
    {
        _egroups.clear();
        _egroups.resize(num_vertices(bg));

        for (auto e : edges_range(g))
        {
            size_t r = b[get_source(e, g)];
            auto& r_elist = _egroups[r];
            _epos[e].first = insert_edge(std::make_tuple(e, true), r_elist,
                                         eweight[e]);

            size_t s = b[get_target(e, g)];
            auto& s_elist = _egroups[s];
            _epos[e].second = insert_edge(std::make_tuple(e, false), s_elist,
                                          eweight[e]);
        }
    }

    void clear()
    {
        _egroups.clear();
    }

    bool empty()
    {
        return _egroups.empty();
    }

    template <class Edge, class EV>
    size_t insert_edge(const Edge& e, EV& elist, size_t)
    {
        elist.push_back(e);
        return elist.size() - 1;
    }

    template <class Edge>
    size_t insert_edge(const Edge& e, DynamicSampler<Edge>& elist,
                       size_t weight)
    {
        return elist.insert(e, weight);
    }

    template <class Edge>
    void remove_edge(size_t pos, vector<Edge>& elist)
    {
        if (get<1>(elist.back()))
            _epos[get<0>(elist.back())].first = pos;
        else
            _epos[get<0>(elist.back())].second = pos;
        if (get<1>(elist[pos]))
            _epos[get<0>(elist[pos])].first = numeric_limits<size_t>::max();
        else
            _epos[get<0>(elist[pos])].second = numeric_limits<size_t>::max();
        elist[pos] = elist.back();
        elist.pop_back();
    }

    template <class Edge>
    void remove_edge(size_t pos, DynamicSampler<Edge>& elist)
    {
        if (get<1>(elist[pos]))
            _epos[get<0>(elist[pos])].first = numeric_limits<size_t>::max();
        else
            _epos[get<0>(elist[pos])].second = numeric_limits<size_t>::max();
        elist.remove(pos);
    }

    template <class Vertex>
    void remove_vertex(Vertex v, Vertex r, Graph& g)
    {
        if (_egroups.empty())
            return;

        typedef Vertex vertex_t;

        auto& elist = _egroups[r];

        // update the half-edge lists
        for (auto e : all_edges_range(v, g))
        {
            vertex_t src = get_source(e, g);
            vertex_t tgt = get_target(e, g);

            bool is_src = (src == v);

            // self-loops will appear twice; we must disambiguate
            if (src == tgt)
            {
                size_t pos = _epos[e].first;
                is_src = (pos < elist.size() && get<0>(elist[pos]) == e);
            }

            size_t pos = (is_src) ? _epos[e].first : _epos[e].second;
            assert(pos < elist.size());
            assert(get<0>(elist[pos]) == e);
            remove_edge(pos, elist);
        }
    }

    template <class Vertex, class Eprop>
    void add_vertex(Vertex v, Vertex s, Eprop& eweight, Graph& g)
    {
        if (_egroups.empty())
            return;

        typedef Vertex vertex_t;

        auto& elist = _egroups[s];

        //update the half-edge lists
        for (auto e : all_edges_range(v, g))
        {
            vertex_t src = get_source(e, g);
            vertex_t tgt = get_target(e, g);

            bool is_src = (src == v);

            // self-loops will appear twice; we must disambiguate
            if (src == tgt)
            {
                size_t pos = _epos[e].first;
                is_src = !(pos < elist.size() && get<0>(elist[pos]) == e);
            }

            typedef typename graph_traits<Graph>::edge_descriptor e_type;
            size_t pos = insert_edge(std::make_tuple(e_type(e), is_src),
                                     elist, size_t(eweight[e]));
            if (is_src)
                _epos[e].first = pos;
            else
                _epos[e].second = pos;
        }
    }

    template <class Vertex, class Eprop>
    void move_vertex(Vertex v, Vertex r, Vertex s, Eprop& eweight,
                     Graph& g)
    {
        if (r == s)
            return;
        remove_egroups(v, r, g);
        add_egroups(v, s, eweight, g);
    }

    template <class Edge, class RNG>
    const auto& sample_edge(const DynamicSampler<Edge>& elist, RNG& rng)
    {
        return get<0>(elist.sample(rng));
    }

    template <class Edge, class RNG>
    const auto& sample_edge(const vector<Edge>& elist, RNG& rng)
    {
        return get<0>(uniform_sample(elist, rng));
    }

    template <class Vertex, class RNG>
    const auto& sample_edge(Vertex r, RNG& rng)
    {
        return sample_edge(_egroups[r], rng);
    }

private:
    typedef typename std::conditional<Weighted::value,
                                      DynamicSampler<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool>>,
                                      vector<std::tuple<typename graph_traits<Graph>::edge_descriptor, bool>>>::type
        sampler_t;
    vector<sampler_t> _egroups;

    typedef typename eprop_map_t<pair<size_t, size_t>>::type epos_t;
    epos_t _epos;
};


// Sample neighbours efficiently
// =============================

template <class Vertex, class Graph, class Eprop, class SMap>
void build_neighbour_sampler(Vertex v, SMap& sampler, Eprop& eweight, Graph& g,
                             bool self_loops=false)
{
    vector<Vertex> neighbours;
    vector<double> probs;
    neighbours.reserve(total_degreeS()(v, g));
    probs.reserve(total_degreeS()(v, g));

    for (auto e : all_edges_range(v, g))
    {
        Vertex u = target(e, g);
        if (is_directed::apply<Graph>::type::value && u == v)
            u = source(e, g);
        if (!self_loops && u == v)
            continue;
        neighbours.push_back(u);
        auto w = eweight[e];
        if (w == 0)
            continue;
        probs.push_back(w);  // Self-loops will be added twice, and hence will
                             // be sampled with probability 2 * eweight[e]
    }

    sampler = Sampler<Vertex, mpl::false_>(neighbours, probs);
};

template <class Vertex, class Graph, class Eprop>
void build_neighbour_sampler(Vertex v, vector<size_t>& sampler, Eprop&, Graph& g,
                             bool self_loops=false)
{
    sampler.clear();
    for (auto e : all_edges_range(v, g))
    {
        Vertex u = target(e, g);
        if (is_directed::apply<Graph>::type::value && u == v)
            u = source(e, g);
        if (!self_loops && u == v)
            continue;
        sampler.push_back(u); // Self-loops will be added twice
    }
};

template <class Graph, class Eprop, class Sampler>
void init_neighbour_sampler(Graph& g, Eprop eweight, Sampler& sampler)
{
    for (auto v : vertices_range(g))
        build_neighbour_sampler(v, sampler[v], eweight, g);
}

template <class Sampler, class RNG>
auto sample_neighbour(Sampler& sampler, RNG& rng)
{
    return sampler.sample(rng);
}

template <class Vertex, class RNG>
auto sample_neighbour(vector<Vertex>& sampler, RNG& rng)
{
    return uniform_sample(sampler, rng);
}


// Sampling marginal probabilities on the edges
template <class Graph, class Vprop, class MEprop>
void collect_edge_marginals(size_t B, Vprop b, MEprop p, Graph& g, Graph&)
{
    for (auto e : edges_range(g))
    {
        auto u = min(source(e, g), target(e, g));
        auto v = max(source(e, g), target(e, g));

        auto r = b[u];
        auto s = b[v];

        auto& pv = p[e];
        if (pv.size() < B * B)
            pv.resize(B * B);
        size_t j = r + B * s;
        pv[j]++;
    }
}

template <class Graph, class Vprop, class VVprop>
void collect_vertex_marginals(Vprop b, VVprop p, Graph& g)
{
    for (auto v : vertices_range(g))
    {
        auto r = b[v];
        auto& pv = p[v];
        if (pv.size() <= size_t(r))
            pv.resize(r + 1);
        pv[r]++;
    }
}

} // graph_tool namespace

#endif //GRAPH_BLOCKMODEL_UTIL_HH
