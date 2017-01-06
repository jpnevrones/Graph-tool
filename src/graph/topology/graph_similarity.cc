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


#include "graph_python_interface.hh"

#include "graph.hh"
#include "graph_filtering.hh"
#include "graph_properties.hh"
#include "graph_selectors.hh"

#include "graph_similarity.hh"

using namespace std;
using namespace boost;
using namespace graph_tool;

size_t similarity(GraphInterface& gi1, GraphInterface& gi2, boost::any label1,
                  boost::any label2)
{
    size_t s = 0;
    gt_dispatch<>()
        (std::bind(get_similarity(), std::placeholders::_1, std::placeholders::_2,
                   std::placeholders::_3, label2, std::ref(s)),
         all_graph_views(),
         all_graph_views(),
         writable_vertex_properties())
        (gi1.get_graph_view(), gi2.get_graph_view(), label1);
    return s;
}

size_t similarity_fast(GraphInterface& gi1, GraphInterface& gi2, boost::any label1,
                       boost::any label2)
{
    size_t s = 0;
    gt_dispatch<boost::mpl::true_>()
        (std::bind(get_similarity_fast(), std::placeholders::_1,
                   std::placeholders::_2, std::placeholders::_3,
                   label2, std::ref(s)),
         all_graph_views(),
         all_graph_views(),
         vertex_integer_properties())
        (gi1.get_graph_view(), gi2.get_graph_view(), label1);
    return s;
}

void export_similarity()
{
    python::def("similarity", &similarity);
    python::def("similarity_fast", &similarity_fast);
};
