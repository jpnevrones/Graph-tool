// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include "graph_filtering.hh"
#include "graph_python_interface.hh"

#include <boost/python.hpp>
#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/undirected_dfs.hpp>

#include "graph.hh"
#include "graph_selectors.hh"
#include "graph_util.hh"

#ifdef HAVE_BOOST_COROUTINE
#include <boost/coroutine/all.hpp>
#endif // HAVE_BOOST_COROUTINE

using namespace std;
using namespace boost;
using namespace graph_tool;


class DFSVisitorWrapper
{
public:
    DFSVisitorWrapper(GraphInterface& gi, python::object vis)
        : _gi(gi), _vis(vis) {}


    template <class Vertex, class Graph>
    void initialize_vertex(Vertex u, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("initialize_vertex")(PythonVertex<Graph>(gp, u));
    }
    template <class Vertex, class Graph>
    void start_vertex(Vertex u, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("start_vertex")(PythonVertex<Graph>(gp, u));
    }
    template <class Vertex, class Graph>
    void discover_vertex(Vertex u, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("discover_vertex")(PythonVertex<Graph>(gp, u));
    }

    template <class Edge, class Graph>
    void examine_edge(Edge e, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("examine_edge")(PythonEdge<Graph>(gp, e));
    }

    template <class Edge, class Graph>
    void tree_edge(Edge e, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("tree_edge")(PythonEdge<Graph>(gp, e));
    }

    template <class Edge, class Graph>
    void back_edge(Edge e, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("back_edge")(PythonEdge<Graph>(gp, e));
    }

    template <class Edge, class Graph>
    void forward_or_cross_edge(Edge e, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("forward_or_cross_edge")(PythonEdge<Graph>(gp, e));
    }

    template <class Vertex, class Graph>
    void finish_vertex(Vertex u, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _vis.attr("finish_vertex")(PythonVertex<Graph>(gp, u));
    }

private:
    GraphInterface& _gi;
    python::object _vis;
};

struct do_dfs
{
    template <class Graph, class VertexIndexMap, class Visitor>
    void operator()(Graph& g, VertexIndexMap vertex_index, size_t s,
                    Visitor vis) const
    {
        typename property_map_type::apply<default_color_type,
                                          VertexIndexMap>::type
            color(vertex_index);
        depth_first_visit(g, vertex(s, g), vis, color);
    }
};


void dfs_search(GraphInterface& g, size_t s, python::object vis)
{
    run_action<graph_tool::all_graph_views,mpl::true_>()
        (g, std::bind(do_dfs(), std::placeholders::_1, g.get_vertex_index(),
                      s, DFSVisitorWrapper(g, vis)))();
}

#ifdef HAVE_BOOST_COROUTINE

class DFSGeneratorVisitor : public dfs_visitor<>
{
public:
    DFSGeneratorVisitor(GraphInterface& gi,
                        coro_t::push_type& yield)
        : _gi(gi), _yield(yield) {}

    template <class Edge, class Graph>
    void tree_edge(const Edge& e, Graph& g)
    {
        auto gp = retrieve_graph_view<Graph>(_gi, g);
        _yield(boost::python::object(PythonEdge<Graph>(gp, e)));
    }

private:
    GraphInterface& _gi;
    coro_t::push_type& _yield;
};

#endif // HAVE_BOOST_COROUTINE


boost::python::object dfs_search_generator(GraphInterface& g, size_t s)
{
#ifdef HAVE_BOOST_COROUTINE
    auto dispatch = [&](auto& yield)
        {
            DFSGeneratorVisitor vis(g, yield);
            run_action<graph_tool::all_graph_views,mpl::true_>()
                (g, std::bind(do_dfs(), std::placeholders::_1,
                              g.get_vertex_index(), s, vis))();
        };
    return boost::python::object(CoroGenerator(dispatch));
#else
    throw GraphException("This functionality is not available because boost::coroutine was not found at compile-time");
#endif
}

void export_dfs()
{
    using namespace boost::python;
    def("dfs_search", &dfs_search);
    def("dfs_search_generator", &dfs_search_generator);
}
