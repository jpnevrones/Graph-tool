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

#ifndef GRAPH_ADJACENCY_HH
#define GRAPH_ADJACENCY_HH

#include <vector>
#include <deque>
#include <utility>
#include <numeric>
#include <iostream>
#include <tuple>
#include <functional>
#include <boost/iterator.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/range/irange.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/properties.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "transform_iterator.hh"

namespace boost
{

// ========================================================================
// Forward declarations
// ========================================================================

template <class Vertex>
class adj_list;

// forward declaration of manipulation functions
template <class Vertex>
std::pair<typename adj_list<Vertex>::vertex_iterator,
          typename adj_list<Vertex>::vertex_iterator>
vertices(const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::edge_iterator,
          typename adj_list<Vertex>::edge_iterator>
edges(const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::edge_descriptor, bool>
edge(Vertex s, Vertex t, const adj_list<Vertex>& g);

template <class Vertex>
size_t out_degree(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
size_t in_degree(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::out_edge_iterator,
          typename adj_list<Vertex>::out_edge_iterator>
out_edges(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::in_edge_iterator,
          typename adj_list<Vertex>::in_edge_iterator>
in_edges(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
adjacent_vertices(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
out_neighbours(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
in_neighbours(Vertex v, const adj_list<Vertex>& g);

template <class Vertex>
size_t num_vertices(const adj_list<Vertex>& g);

template <class Vertex>
size_t num_edges(const adj_list<Vertex>& g);

template <class Vertex>
Vertex add_vertex(adj_list<Vertex>& g);

template <class Vertex>
void clear_vertex(Vertex v, adj_list<Vertex>& g);

template <class Vertex>
void remove_vertex(Vertex v, adj_list<Vertex>& g);

template <class Vertex>
void remove_vertex_fast(Vertex v, adj_list<Vertex>& g);

template <class Vertex>
std::pair<typename adj_list<Vertex>::edge_descriptor, bool>
add_edge(Vertex s, Vertex t, adj_list<Vertex>& g);

template <class Vertex>
void remove_edge(Vertex s, Vertex t, adj_list<Vertex>& g);

template <class Vertex>
void remove_edge(const typename adj_list<Vertex>::edge_descriptor& e,
                 adj_list<Vertex>& g);

// ========================================================================
// adj_list<Vertex>
// ========================================================================
//
// adj_list is a very simple adjacency list implementation for bidirectional
// graphs based on std::vector, meant to be reasonably efficient both
// memory-wise and computationally. It maintains a list of in and out-edges for
// each vertex, and each edge has a built-in index (which is replicated in both
// lists). For each edge, a total of 4 integers is necessary: the source and
// target vertices, in the in_edges and out_edges lists, respectively, and the
// (same) edge index in both lists. The integer type is given by the Vertex
// template parameter. It achieves about half as much memory as
// boost::adjacency_list with an edge index property map and the same integer
// type.

// The complexity guarantees and iterator invalidation rules are the same as
// boost::adjacency_list with vector storage selectors for both vertex and edge
// lists.

namespace detail
{
template <class Vertex>
struct adj_edge_descriptor
{
    adj_edge_descriptor()
        : s(std::numeric_limits<Vertex>::max()),
          t(std::numeric_limits<Vertex>::max()),
          idx(std::numeric_limits<Vertex>::max()), inv(false) {};
    adj_edge_descriptor(Vertex s, Vertex t, Vertex idx, bool inv)
        : s(s), t(t), idx(idx), inv(inv) {}

    bool operator==(const adj_edge_descriptor& other) const
    {
        return idx == other.idx;
    }
    bool operator!=(const adj_edge_descriptor& other) const
    {
        return idx != other.idx;
    }
    bool operator<(const adj_edge_descriptor& other) const
    {
        return idx < other.idx;
    }

    Vertex s, t, idx;
    bool inv;
};
} // namespace detail

template <class Vertex = size_t>
class adj_list
{
public:
    struct graph_tag {};
    typedef Vertex vertex_t;

    typedef detail::adj_edge_descriptor<Vertex> edge_descriptor;

    typedef std::vector<std::pair<vertex_t, vertex_t> > edge_list_t;
    typedef std::vector<edge_list_t> vertex_list_t;
    typedef typename integer_range<Vertex>::iterator vertex_iterator;

    adj_list(): _n_edges(0), _edge_index_range(0), _keep_epos(false) {}

    struct get_vertex
    {
        get_vertex() {}
        typedef Vertex result_type;
        __attribute__((always_inline))
        Vertex operator()(const std::pair<vertex_t, vertex_t>& v) const
        { return v.first; }
    };

    typedef transform_random_access_iterator<get_vertex,
                                             typename edge_list_t::const_iterator>
        adjacency_iterator;

    typedef adjacency_iterator in_adjacency_iterator;

    template <class Deference>
    struct base_edge_iterator:
        public boost::iterator_facade<base_edge_iterator<Deference>,
                                      edge_descriptor,
                                      std::random_access_iterator_tag,
                                      edge_descriptor>
    {
        base_edge_iterator() {}
        base_edge_iterator(vertex_t v, typename edge_list_t::const_iterator&& iter)
            : _v(v), _iter(std::forward<typename edge_list_t::const_iterator>(iter))
        {}

    private:
        friend class boost::iterator_core_access;
        void increment() { ++_iter; }
        void decrement() { --_iter; }
        template <class Distance>
        void advance(Distance n) { _iter += n; }
        auto distance_to(base_edge_iterator const& other) const
        {
            return other._iter - _iter;
        }

        bool equal(base_edge_iterator const& other) const
        {
            return _iter == other._iter;
        }

        edge_descriptor dereference() const
        {
            return Deference::def(_v, *_iter);
        }

        vertex_t _v;
        typename edge_list_t::const_iterator _iter;
    };

    struct make_out_edge
    {
        static edge_descriptor def(vertex_t src,
                                   const std::pair<vertex_t, vertex_t>& v)
        { return edge_descriptor(src, v.first, v.second, false); }
    };

    struct make_in_edge
    {
        static edge_descriptor def(vertex_t tgt,
                                   const std::pair<vertex_t, vertex_t>& v)
        { return edge_descriptor(v.first, tgt, v.second, false); }
    };

    typedef base_edge_iterator<make_out_edge> out_edge_iterator;
    typedef base_edge_iterator<make_in_edge> in_edge_iterator;

    class edge_iterator:
        public boost::iterator_facade<edge_iterator,
                                      edge_descriptor,
                                      boost::forward_traversal_tag,
                                      edge_descriptor>
    {
    public:
        edge_iterator() {}
        explicit edge_iterator(const typename vertex_list_t::const_iterator& vi_begin,
                               const typename vertex_list_t::const_iterator& vi_end,
                               const typename vertex_list_t::const_iterator& vi,
                               const typename edge_list_t::const_iterator& ei)
            : _vi_begin(vi_begin), _vi_end(vi_end), _vi(vi), _ei(ei)
        {
            // move position to first edge
            skip();
        }

    private:
        friend class boost::iterator_core_access;

        void skip()
        {
            //skip empty vertices
            while (_vi != _vi_end && _ei == _vi->end())
            {
                ++_vi;
                if (_vi != _vi_end)
                    _ei = _vi->begin();
            }
        }

        void increment()
        {
            ++_ei;
            skip();
        }

        bool equal(edge_iterator const& other) const
        {
            if (_vi_begin == _vi_end)
                return _vi == other._vi;
            return _vi == other._vi && _ei == other._ei;
        }

        edge_descriptor dereference() const
        {
            return edge_descriptor(vertex_t(_vi - _vi_begin),
                                   _ei->first, _ei->second, false);
        }

        typename vertex_list_t::const_iterator _vi_begin;
        typename vertex_list_t::const_iterator _vi_end;
        typename vertex_list_t::const_iterator _vi;
        typename edge_list_t::const_iterator _ei;
    };

    void reindex_edges()
    {
        _free_indexes.clear();
        _edge_index_range = 0;
        _in_edges.clear();
        _in_edges.resize(_out_edges.size());
        for (size_t i = 0; i < _out_edges.size(); ++i)
        {
            auto& oes = _out_edges[i];
            for (size_t j = 0; j < oes.size(); ++j)
            {
                auto& oe = oes[j];
                Vertex v = oe.first;
                oe.second = _edge_index_range;
                _in_edges[v].emplace_back(i, _edge_index_range);
                _edge_index_range++;
            }
        }

        if (_keep_epos)
            rebuild_epos();
    }

    void set_keep_epos(bool keep)
    {
        if (keep)
        {
            if (!_keep_epos)
                rebuild_epos();
        }
        else
        {
            _epos.clear();
        }
        _keep_epos = keep;
    }

    bool get_keep_epos()
    {
        return _keep_epos;
    }

    size_t get_edge_index_range() const { return _edge_index_range; }

    static Vertex null_vertex() { return std::numeric_limits<Vertex>::max(); }

    void shrink_to_fit()
    {
        _in_edges.shrink_to_fit();
        _out_edges.shrink_to_fit();
        std::for_each(_in_edges.begin(), _in_edges.end(),
                      [](auto &es){es.shrink_to_fit();});
        std::for_each(_out_edges.begin(), _out_edges.end(),
                      [](auto &es){es.shrink_to_fit();});
        auto erange = boost::edges(*this);
        auto iter = std::max_element(erange.first, erange.second,
                                     [](const auto &a, const auto& b) -> bool
                                     {return a.idx < b.idx;});
        if (iter == erange.second)
            _edge_index_range = 0;
        else
            _edge_index_range = iter->idx + 1;
        auto iter_idx = std::remove_if(_free_indexes.begin(),
                                       _free_indexes.end(),
                                       [&](auto idx) -> bool
                                       {return idx > _edge_index_range;});
        _free_indexes.erase(iter_idx, _free_indexes.end());
        _free_indexes.shrink_to_fit();
        if (_keep_epos)
            _epos.resize(_edge_index_range);
        _epos.shrink_to_fit();
    }

private:
    vertex_list_t _out_edges;
    vertex_list_t _in_edges;
    size_t _n_edges;
    size_t _edge_index_range;
    std::deque<size_t> _free_indexes; // indexes of deleted edges to be used up
                                      // for new edges to avoid very large
                                      // indexes, and unnecessary property map
                                      // memory use
    bool _keep_epos;
    std::vector<std::pair<int32_t, int32_t>> _epos;

    void rebuild_epos()
    {
        _epos.resize(_edge_index_range);
        for (size_t i = 0; i < _out_edges.size(); ++i)
        {
            auto& oes = _out_edges[i];
            for (size_t j = 0; j < oes.size(); ++j)
            {
                size_t idx = oes[j].second;
                _epos[idx].first = j;
            }

            auto& ies = _in_edges[i];
            for (size_t j = 0; j < ies.size(); ++j)
            {
                size_t idx = _in_edges[i][j].second;
                _epos[idx].second = j;
            }
        }
    }

    // manipulation functions
    friend std::pair<vertex_iterator, vertex_iterator>
    vertices<>(const adj_list<Vertex>& g);

    friend std::pair<edge_iterator, edge_iterator>
    edges<>(const adj_list<Vertex>& g);

    friend std::pair<edge_descriptor, bool>
    edge<>(Vertex s, Vertex t, const adj_list<Vertex>& g);

    friend size_t out_degree<>(Vertex v, const adj_list<Vertex>& g);

    friend size_t in_degree<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<out_edge_iterator, out_edge_iterator>
    out_edges<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<in_edge_iterator, in_edge_iterator>
    in_edges<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<adjacency_iterator, adjacency_iterator>
    adjacent_vertices<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<adjacency_iterator, adjacency_iterator>
    out_neighbours<>(Vertex v, const adj_list<Vertex>& g);

    friend std::pair<adjacency_iterator, adjacency_iterator>
    in_neighbours<>(Vertex v, const adj_list<Vertex>& g);

    friend size_t num_vertices<>(const adj_list<Vertex>& g);

    friend size_t num_edges<>(const adj_list<Vertex>& g);

    friend Vertex add_vertex<>(adj_list<Vertex>& g);

    friend void clear_vertex<>(Vertex v, adj_list<Vertex>& g);

    friend void remove_vertex<>(Vertex v, adj_list<Vertex>& g);

    friend void remove_vertex_fast<>(Vertex v, adj_list<Vertex>& g);

    friend std::pair<edge_descriptor, bool>
    add_edge<>(Vertex s, Vertex t, adj_list<Vertex>& g);

    friend void remove_edge<>(Vertex s, Vertex t, adj_list<Vertex>& g);

    friend void remove_edge<>(const edge_descriptor& e, adj_list<Vertex>& g);
};

//========================================================================
// Graph traits and BGL scaffolding
//========================================================================

struct adj_list_traversal_tag
    : public vertex_list_graph_tag,
      public edge_list_graph_tag,
      public adjacency_graph_tag,
      public bidirectional_graph_tag,
      public adjacency_matrix_tag {};

template <class Vertex>
struct graph_traits<adj_list<Vertex> >
{
    typedef Vertex vertex_descriptor;
    typedef typename adj_list<Vertex>::edge_descriptor edge_descriptor;
    typedef typename adj_list<Vertex>::edge_iterator edge_iterator;
    typedef typename adj_list<Vertex>::adjacency_iterator adjacency_iterator;

    typedef typename adj_list<Vertex>::out_edge_iterator out_edge_iterator;
    typedef typename adj_list<Vertex>::in_edge_iterator in_edge_iterator;

    typedef typename adj_list<Vertex>::vertex_iterator vertex_iterator;

    typedef bidirectional_tag directed_category;
    typedef allow_parallel_edge_tag edge_parallel_category;
    typedef adj_list_traversal_tag traversal_category;

    typedef Vertex vertices_size_type;
    typedef Vertex edges_size_type;
    typedef size_t degree_size_type;

    static Vertex null_vertex() { return adj_list<Vertex>::null_vertex(); }

private:
    BOOST_STATIC_ASSERT((is_convertible<typename std::iterator_traits<out_edge_iterator>::iterator_category,
                                        std::random_access_iterator_tag>::value));
    BOOST_STATIC_ASSERT((is_convertible<typename std::iterator_traits<in_edge_iterator>::iterator_category,
                                        std::random_access_iterator_tag>::value));
    BOOST_STATIC_ASSERT((is_convertible<typename std::iterator_traits<adjacency_iterator>::iterator_category,
                                        std::random_access_iterator_tag>::value));
};

template <class Vertex>
struct graph_traits<const adj_list<Vertex> >
    : public graph_traits<adj_list<Vertex> >
{
};


template <class Vertex>
struct edge_property_type<adj_list<Vertex> >
{
    typedef void type;
};

template <class Vertex>
struct vertex_property_type<adj_list<Vertex> >
{
    typedef void type;
};

template <class Vertex>
struct graph_property_type<adj_list<Vertex> >
{
    typedef void type;
};

//========================================================================
// Graph access and manipulation functions
//========================================================================

template <class Vertex>
inline __attribute__((always_inline))
std::pair<typename adj_list<Vertex>::vertex_iterator,
          typename adj_list<Vertex>::vertex_iterator>
vertices(const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::vertex_iterator vi_t;
    return std::make_pair(vi_t(0), vi_t(g._out_edges.size()));
}


template <class Vertex>
inline
std::pair<typename adj_list<Vertex>::edge_iterator,
          typename adj_list<Vertex>::edge_iterator>
edges(const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::edge_list_t::const_iterator ei_t;
    typedef typename adj_list<Vertex>::vertex_list_t::const_iterator vi_t;
    ei_t ei_begin, ei_end;
    vi_t last_vi;
    if (g._out_edges.empty())
    {
        last_vi = g._out_edges.end();
    }
    else
    {
        ei_begin = g._out_edges[0].begin();
        last_vi = g._out_edges.end() - 1;
        ei_end = last_vi->end();
    }
    typename adj_list<Vertex>::edge_iterator ebegin(g._out_edges.begin(),
                                                    g._out_edges.end(),
                                                    g._out_edges.begin(),
                                                    ei_begin);
    typename adj_list<Vertex>::edge_iterator eend(g._out_edges.begin(),
                                                  g._out_edges.end(),
                                                  last_vi,
                                                  ei_end);
    return std::make_pair(ebegin, eend);
}

template <class Vertex>
inline __attribute__((always_inline))
Vertex vertex(size_t i, const adj_list<Vertex>&)
{
    return i;
}

template <class Vertex>
inline
std::pair<typename adj_list<Vertex>::edge_descriptor, bool>
edge(Vertex s, Vertex t, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::edge_descriptor edge_descriptor;
    const auto& oes = g._out_edges[s];
    auto iter = std::find_if(oes.begin(), oes.end(),
                             [&](const auto& e) -> bool {return e.first == t;});
    if (iter != oes.end())
        return std::make_pair(edge_descriptor(s, t, iter->second, false),
                              true);
    Vertex v = graph_traits<adj_list<Vertex> >::null_vertex();
    return std::make_pair(edge_descriptor(v, v, v, false), false);
}

template <class Vertex>
inline __attribute__((always_inline))
size_t out_degree(Vertex v, const adj_list<Vertex>& g)
{
    return g._out_edges[v].size();
}

template <class Vertex>
inline __attribute__((always_inline))
size_t in_degree(Vertex v, const adj_list<Vertex>& g)
{
    return g._in_edges[v].size();
}

template <class Vertex>
inline __attribute__((always_inline))
size_t degree(Vertex v, const adj_list<Vertex>& g)
{
    return in_degree(v, g) + out_degree(v, g);
}

template <class Vertex>
inline __attribute__((always_inline))
std::pair<typename adj_list<Vertex>::out_edge_iterator,
          typename adj_list<Vertex>::out_edge_iterator>
out_edges(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::out_edge_iterator ei_t;
    auto& edges = g._out_edges[v];
    return std::make_pair(ei_t(v, edges.begin()),
                          ei_t(v, edges.end()));
}

template <class Vertex>
inline  __attribute__((always_inline))
std::pair<typename adj_list<Vertex>::in_edge_iterator,
          typename adj_list<Vertex>::in_edge_iterator>
in_edges(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::in_edge_iterator ei_t;
    auto& edges = g._in_edges[v];
    return std::make_pair(ei_t(v, edges.begin()),
                          ei_t(v, edges.end()));
}

template <class Vertex>
inline __attribute__((always_inline))
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
out_neighbours(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::adjacency_iterator ai_t;
    auto& edges = g._out_edges[v];
    return std::make_pair(ai_t(edges.begin()),
                          ai_t(edges.end()));
}

template <class Vertex>
inline __attribute__((always_inline))
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
in_neighbours(Vertex v, const adj_list<Vertex>& g)
{
    typedef typename adj_list<Vertex>::adjacency_iterator ai_t;
    auto& edges = g._in_edges[v];
    return std::make_pair(ai_t(edges.begin()),
                          ai_t(edges.end()));
}

template <class Vertex>
inline __attribute__((always_inline))
std::pair<typename adj_list<Vertex>::adjacency_iterator,
          typename adj_list<Vertex>::adjacency_iterator>
adjacent_vertices(Vertex v, const adj_list<Vertex>& g)
{
    return out_neighbours(v, g);
}


template <class Vertex>
inline __attribute__((always_inline))
size_t num_vertices(const adj_list<Vertex>& g)
{
    return g._out_edges.size();
}

template <class Vertex>
inline __attribute__((always_inline))
size_t num_edges(const adj_list<Vertex>& g)
{
    return g._n_edges;
}

template <class Vertex>
inline __attribute__((always_inline))
Vertex add_vertex(adj_list<Vertex>& g)
{
    g._out_edges.emplace_back();
    g._in_edges.emplace_back();
    return g._out_edges.size() - 1;
}

template <class Vertex>
inline void clear_vertex(Vertex v, adj_list<Vertex>& g)
{
    if (!g._keep_epos)
    {
        auto remove_es = [&] (auto& out_edges, auto& in_edges)
            {
                auto& oes = out_edges[v];
                for (const auto& oe : oes)
                {
                    Vertex t = oe.first;
                    auto& ies = in_edges[t];
                    auto iter =
                        std::remove_if(ies.begin(), ies.end(),
                                       [&](const auto& ei) -> bool
                                       {
                                           if (ei.first == v)
                                           {
                                               g._free_indexes.push_back(ei.second);
                                               return true;
                                           }
                                           return false;
                                       });
                    ies.erase(iter, ies.end());
                }
                g._n_edges -= oes.size();
                oes.clear();
            };
        remove_es(g._out_edges, g._in_edges);
        remove_es(g._in_edges, g._out_edges);
    }
    else
    {
        auto remove_es = [&] (auto& out_edges, auto& in_edges,
                              const auto& get_pos)
        {
            auto& oes = out_edges[v];
            for (const auto& ei : oes)
            {
                Vertex t = ei.first;
                size_t idx = ei.second;
                auto& ies = in_edges[t];
                auto& back = ies.back();
                auto& pos = get_pos(idx);
                auto& bpos = get_pos(back.second);
                bpos = pos;
                ies[pos] = back;
                ies.pop_back();
                g._free_indexes.push_back(idx);
            }
            g._n_edges -= oes.size();
            oes.clear();
        };
        remove_es(g._out_edges, g._in_edges,
                  [&](size_t idx) -> auto& {return g._epos[idx].second;});
        remove_es(g._in_edges, g._out_edges,
                  [&](size_t idx) -> auto& {return g._epos[idx].first;});
    }
}

// O(V + E)
template <class Vertex>
inline void remove_vertex(Vertex v, adj_list<Vertex>& g)
{
    clear_vertex(v, g);
    g._out_edges.erase(g._out_edges.begin() + v);
    g._in_edges.erase(g._in_edges.begin() + v);

    auto shift_es = [&](auto& edges, int i)
    {
        for (auto& e : edges[i])
        {
            if (e.first > v)
                e.first--;
        }
    };

    size_t N = g._out_edges.size();
    #pragma omp parallel for schedule(runtime) if (N > 100)
    for (size_t i = 0; i < N; ++i)
    {
        shift_es(g._out_edges, i);
        shift_es(g._in_edges, i);
    }
}

// O(k + k_last)
template <class Vertex>
inline void remove_vertex_fast(Vertex v, adj_list<Vertex>& g)
{
    Vertex back = g._out_edges.size() - 1;

    if (v < back)
    {
        clear_vertex(v, g);
        g._out_edges[v].swap(g._out_edges[back]);
        g._in_edges[v].swap(g._in_edges[back]);
        g._out_edges.pop_back();
        g._in_edges.pop_back();

        auto rename_v = [&] (auto& out_edges, auto& in_edges,
                             const auto& get_pos)
            {
                for (auto& eu : out_edges[v])
                {
                    Vertex u = eu.first;
                    if (u == back)
                    {
                        eu.first = v;
                    }
                    else
                    {
                        if (!g._keep_epos)
                        {
                            for (auto& e : in_edges[u])
                            {
                                if (e.first == back)
                                    e.first = v;
                            }
                        }
                        else
                        {
                            size_t idx = eu.second;
                            auto pos = get_pos(idx);
                            in_edges[u][pos].first = v;
                        }
                    }
                }
            };

        rename_v(g._out_edges, g._in_edges,
                 [&](size_t idx) -> auto {return g._epos[idx].second;});
        rename_v(g._in_edges, g._out_edges,
                 [&](size_t idx) -> auto {return g._epos[idx].first;});
    }
    else
    {
        clear_vertex(v, g);
        g._out_edges.pop_back();
        g._in_edges.pop_back();
    }
}

template <class Vertex>
inline
typename std::pair<typename adj_list<Vertex>::edge_descriptor, bool>
add_edge(Vertex s, Vertex t, adj_list<Vertex>& g)
{
    Vertex idx;
    if (g._free_indexes.empty())
    {
        idx = g._edge_index_range++;
    }
    else
    {
        idx = g._free_indexes.front();
        g._free_indexes.pop_front();
    }

    auto& oes = g._out_edges[s];
    auto& ies = g._in_edges[t];
    oes.emplace_back(t, idx);
    ies.emplace_back(s, idx);
    g._n_edges++;

    if (g._keep_epos)
    {
        if (idx >= g._epos.size())
            g._epos.resize(idx + 1);
        auto& ei = g._epos[idx];
        ei.first = oes.size() - 1;
        ei.second = ies.size() - 1;
    }

    typedef typename adj_list<Vertex>::edge_descriptor edge_descriptor;
    return std::make_pair(edge_descriptor(s, t, idx, false), true);
}

template <class Vertex>
inline void remove_edge(Vertex s, Vertex t,
                        adj_list<Vertex>& g)
{
    if (!g._keep_epos)
    {
        auto& oes = g._out_edges[s];
        auto iter_o = std::find_if(oes.begin(), oes.end(),
                                   [&] (const auto& ei) -> bool
                                   {return t == ei.first;});
        if (iter_o != oes.end())
        {
            g._free_indexes.push_back(iter_o->second);
            oes.erase(iter_o);
            g._n_edges--;
        }

        auto& ies = g._in_edges[t];
        auto iter_i = std::find_if(ies.begin(), ies.end(),
                                   [&] (const auto& ei) -> bool
                                   {return s == ei.first;});
        if (iter_i != ies.end())
        {
            ies.erase(iter_i);
        }
    }
    else
    {
        remove_edge(edge(s, t, g).first, g);
    }
}

template <class Vertex>
inline void remove_edge(const typename adj_list<Vertex>::edge_descriptor& e,
                        adj_list<Vertex>& g)
{
    auto& s = e.s;
    auto& t = e.t;
    auto& idx = e.idx;
    auto& oes = g._out_edges[s];
    auto& ies = g._in_edges[t];

    bool found = false;
    if (!g._keep_epos) // O(k_s + k_t)
    {
        auto remove_e = [&] (auto& elist, auto v)
            {
                auto iter = std::find_if(elist.begin(), elist.end(),
                                         [&] (const auto& ei) -> bool
                                         {return v == ei.first && idx == ei.second;});
                if (iter != elist.end())
                {
                    elist.erase(iter);
                    found = true;
                }
            };

        remove_e(oes, t);
        remove_e(ies, s);
    }
    else // O(1)
    {
        if (idx < g._epos.size())
        {
            auto remove_e = [&] (auto& elist, const auto& get_pos)
            {
                const auto& back = elist.back();
                auto pindex = get_pos(idx);
                get_pos(back.second) = pindex;
                elist[pindex] = back;
                elist.pop_back();
            };

            remove_e(oes, [&](size_t idx) -> auto& {return g._epos[idx].first;});
            remove_e(ies, [&](size_t idx) -> auto& {return g._epos[idx].second;});

            found = true;
        }
    }

    if (found)
    {
        g._free_indexes.push_back(idx);
        g._n_edges--;
    }
}


template <class Vertex>
inline
Vertex source(const typename adj_list<Vertex>::edge_descriptor& e,
              const adj_list<Vertex>&)
{
    return e.s;
}

template <class Vertex>
inline
Vertex target(const typename adj_list<Vertex>::edge_descriptor& e,
              const adj_list<Vertex>&)
{
    return e.t;
}

//========================================================================
// Vertex and edge index property maps
//========================================================================

template <class Vertex>
struct property_map<adj_list<Vertex>, vertex_index_t>
{
    typedef identity_property_map type;
    typedef type const_type;
};

template <class Vertex>
struct property_map<const adj_list<Vertex>, vertex_index_t>
{
    typedef identity_property_map type;
    typedef type const_type;
};

template <class Vertex>
inline identity_property_map
get(vertex_index_t, adj_list<Vertex>&)
{
    return identity_property_map();
}

template <class Vertex>
inline identity_property_map
get(vertex_index_t, const adj_list<Vertex>&)
{
    return identity_property_map();
}

template<class Vertex>
class adj_edge_index_property_map:
    public put_get_helper<Vertex, adj_edge_index_property_map<Vertex> >
{
public:
    typedef typename adj_list<Vertex>::edge_descriptor key_type;
    typedef Vertex reference;
    typedef Vertex value_type;
    typedef boost::readable_property_map_tag category;

    reference operator[](const key_type& k) const {return k.idx;}
};

template <class Vertex>
struct property_map<adj_list<Vertex>, edge_index_t>
{
    typedef adj_edge_index_property_map<Vertex> type;
    typedef type const_type;

};

template <class Vertex>
inline adj_edge_index_property_map<Vertex>
get(edge_index_t, const adj_list<Vertex>&)
{
    return adj_edge_index_property_map<Vertex>();
}

} // namespace boost

// hashing of edge descriptors

namespace std
{

template <class Vertex>
struct hash<boost::detail::adj_edge_descriptor<Vertex>>
{
    template <class Edge>
    std::size_t operator()(Edge const& e) const
    {
        return _h(e.idx);
    }
    std::hash<Vertex> _h;
};

} // namespace std


#endif //GRAPH_ADJACENCY_HH
