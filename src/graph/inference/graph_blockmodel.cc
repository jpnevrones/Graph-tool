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

#include "graph_tool.hh"
#include "random.hh"

#include <boost/python.hpp>

#include "graph_blockmodel_util.hh"
#include "graph_blockmodel.hh"
#include "graph_state.hh"

using namespace boost;
using namespace graph_tool;

namespace graph_tool
{
vmap_t::unchecked_t uncheck(boost::any& amap, vmap_t::unchecked_t*)
{
    return any_cast<vmap_t&>(amap).get_unchecked();
}
emap_t::unchecked_t uncheck(boost::any& amap, emap_t::unchecked_t*)
{
    return any_cast<emap_t&>(amap).get_unchecked();
}
}

GEN_DISPATCH(block_state, BlockState, BLOCK_STATE_params)

python::object make_block_state(boost::python::object ostate,
                                rng_t& rng);

degs_map_t get_block_degs(GraphInterface& gi, boost::any ab, boost::any aweight)
{
    degs_map_t degs;
    vmap_t b = boost::any_cast<vmap_t>(ab);
    run_action<>()(gi,
                   [&](auto& g, auto& eweight)
                   {
                       std::vector<gt_hash_map<std::tuple<size_t, size_t>,
                                               size_t>> hist;
                       for (auto v : vertices_range(g))
                       {
                           size_t r = b[v];
                           if (r >= hist.size())
                               hist.resize(r + 1);
                           size_t kin = in_degreeS()(v, g, eweight);
                           size_t kout = out_degreeS()(v, g, eweight);
                           hist[r][std::make_tuple(kin, kout)]++;
                       }

                       for (size_t r = 0; r < hist.size(); ++r)
                       {
                           auto& deg = degs[r];
                           for (auto& kn : hist[r])
                               deg.emplace_back(get<0>(kn.first),
                                                get<1>(kn.first),
                                                kn.second);
                       }
                   },
                   eweight_tr())(aweight);
    return degs;
}

degs_map_t get_weighted_block_degs(GraphInterface& gi, degs_map_t& degs,
                                   boost::any ab)
{
    degs_map_t ndegs;
    vmap_t b = boost::any_cast<vmap_t>(ab);
    run_action<>()(gi,
                   [&](auto& g)
                   {
                       std::vector<gt_hash_map<std::tuple<size_t, size_t>,
                                               size_t>> hist;
                       for (auto v : vertices_range(g))
                       {
                           size_t r = b[v];
                           if (r >= hist.size())
                               hist.resize(r + 1);
                           auto& h = hist[r];
                           auto& ks = degs[v];
                           for (auto& k : ks)
                               h[std::make_tuple(get<0>(k), get<1>(k))] += get<2>(k);
                       }

                       for (size_t r = 0; r < hist.size(); ++r)
                       {
                           auto& deg = ndegs[r];
                           for (auto& kn : hist[r])
                               deg.emplace_back(get<0>(kn.first),
                                                get<1>(kn.first),
                                                kn.second);
                       }
                   })();
    return ndegs;
}


template <class Prop>
boost::any get_any(Prop& p)
{
    return any(p);
}

void print_degs(degs_map_t& degs, size_t B)
{
    for (size_t r = 0; r < B; ++r)
    {
        cout << r << ":: ";
        auto& ks = degs[r];
        for (auto& k : ks)
        {
            cout << "(" << get<0>(k) << ", " << get<1>(k) << "): "
                 << get<2>(k) << "  ";
        }
        cout << endl;
    }
}

degs_map_t copy_degs(degs_map_t& degs)
{
    return degs.copy();
}

simple_degs_t copy_simple_degs(simple_degs_t& degs)
{
    return degs;
}

void export_blockmodel_state()
{
    using namespace boost::python;

    block_state::dispatch
        ([&](auto* s)
         {
             typedef typename std::remove_reference<decltype(*s)>::type state_t;
             void (state_t::*remove_vertex)(size_t) =
                 &state_t::remove_vertex;
             void (state_t::*add_vertex)(size_t, size_t) =
                 &state_t::add_vertex;
             void (state_t::*remove_vertices)(python::object) =
                 &state_t::remove_vertices;
             void (state_t::*add_vertices)(python::object, python::object) =
                 &state_t::add_vertices;
             void (state_t::*move_vertices)(python::object, python::object) =
                 &state_t::move_vertices;
             double (state_t::*virtual_move)(size_t, size_t, size_t,
                                             entropy_args_t) =
                 &state_t::virtual_move;
             size_t (state_t::*sample_block)(size_t, double, rng_t&) =
                 &state_t::sample_block;
             size_t (state_t::*random_neighbour)(size_t, rng_t&) =
                 &state_t::random_neighbour;
             double (state_t::*get_move_prob)(size_t, size_t, size_t, double,
                                              bool) =
                 &state_t::get_move_prob;
             void (state_t::*merge_vertices)(size_t, size_t) =
                 &state_t::merge_vertices;
             void (state_t::*set_partition)(boost::any&) =
                 &state_t::set_partition;

             class_<state_t> c(name_demangle(typeid(state_t).name()).c_str(),
                               no_init);
             c.def("remove_vertex", remove_vertex)
                 .def("add_vertex", add_vertex)
                 .def("remove_vertices", remove_vertices)
                 .def("add_vertices", add_vertices)
                 .def("move_vertex", &state_t::move_vertex)
                 .def("move_vertices", move_vertices)
                 .def("set_partition", set_partition)
                 .def("virtual_move", virtual_move)
                 .def("merge_vertices", merge_vertices)
                 .def("sample_block", sample_block)
                 .def("sample_neighbour", random_neighbour)
                 .def("entropy", &state_t::entropy)
                 .def("get_partition_dl", &state_t::get_partition_dl)
                 .def("get_deg_dl", &state_t::get_deg_dl)
                 .def("get_move_prob", get_move_prob)
                 .def("enable_partition_stats",
                      &state_t::enable_partition_stats)
                 .def("disable_partition_stats",
                      &state_t::disable_partition_stats)
                 .def("is_partition_stats_enabled",
                      &state_t::is_partition_stats_enabled)
                 .def("couple_state",
                      &state_t::couple_state)
                 .def("decouple_state",
                      &state_t::decouple_state)
                 .def("clear_egroups",
                      &state_t::clear_egroups)
                 .def("rebuild_neighbour_sampler",
                      &state_t::rebuild_neighbour_sampler)
                 .def("sync_emat",
                      &state_t::sync_emat);
         });

    class_<vcmap_t>("unity_vprop_t").def("_get_any", &get_any<vcmap_t>);
    class_<ecmap_t>("unity_eprop_t").def("_get_any", &get_any<ecmap_t>);

    class_<entropy_args_t>("entropy_args")
        .def_readwrite("exact", &entropy_args_t::exact)
        .def_readwrite("dense", &entropy_args_t::dense)
        .def_readwrite("multigraph", &entropy_args_t::multigraph)
        .def_readwrite("adjacency", &entropy_args_t::adjacency)
        .def_readwrite("partition_dl", &entropy_args_t::partition_dl)
        .def_readwrite("degree_dl", &entropy_args_t::degree_dl)
        .def_readwrite("degree_dl_kind", &entropy_args_t::degree_dl_kind)
        .def_readwrite("edges_dl", &entropy_args_t::edges_dl);

    enum_<deg_dl_kind>("deg_dl_kind")
        .value("ent", deg_dl_kind::ENT)
        .value("uniform", deg_dl_kind::UNIFORM)
        .value("dist", deg_dl_kind::DIST);

    enum_<weight_type>("rec_type")
        .value("none", weight_type::NONE)
        .value("positive", weight_type::POSITIVE)
        .value("signed", weight_type::SIGNED)
        .value("delta_t", weight_type::DELTA_T);

    def("make_block_state", &make_block_state);

    def("get_block_degs", &get_block_degs);
    def("get_weighted_block_degs", &get_weighted_block_degs);
    class_<degs_map_t>("degs_map_t")
        .def("print", &print_degs)
        .def("copy", &copy_degs);
    class_<simple_degs_t>("simple_degs_t")
        .def("copy", &copy_simple_degs);

    def("init_q_cache", init_q_cache);
    def("log_q", log_q<size_t>);
    def("q_rec", q_rec);
    def("q_rec_memo", q_rec_memo);
    def("log_q_approx", log_q_approx);
    def("log_q_approx_big", log_q_approx_big);
    def("log_q_approx_small", log_q_approx_small);
}
