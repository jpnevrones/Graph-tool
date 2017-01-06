#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
# graph_tool -- a general graph manipulation python module
#
# Copyright (C) 2006-2016 Tiago de Paula Peixoto <tiago@skewed.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import division, absolute_import, print_function
import sys
if sys.version_info < (3,):
    range = xrange

from .. import _degree, _prop, Graph, GraphView, libcore, _get_rng, PropertyMap, \
    conv_pickle_state, Vector_size_t, Vector_double, group_vector_property
from .. generation import condensation_graph
from .. stats import label_self_loops
from .. spectral import adjacency
import random
from numpy import *
import numpy
import copy
import collections

from . util import *

from .. dl_import import dl_import
dl_import("from . import libgraph_tool_inference as libinference")

from . libgraph_tool_inference import PartitionHist, BlockPairHist

__test__ = False

def set_test(test):
    global __test__
    __test__ = test

def _bm_test():
    global __test__
    return __test__

def get_block_graph(g, B, b, vcount=None, ecount=None, rec=None, drec=None):
    if isinstance(ecount, libinference.unity_eprop_t):
        ecount = None
    if isinstance(vcount, libinference.unity_vprop_t):
        vcount = None
    avprops = []
    if vcount is not None:
        avprops.append(vcount)
    aeprops = []
    if ecount is not None:
        aeprops.append(ecount)
    if rec is not None:
        aeprops.append(rec)
    if drec is not None:
        aeprops.append(drec)
    cg, br, vc, ec, av, ae = condensation_graph(g, b,
                                                avprops=avprops,
                                                aeprops=aeprops,
                                                self_loops=True)
    if vcount is not None:
        vcount = av[0]
        del av[0]
    else:
        vcount = vc
    cg.vp.count = vcount

    if ecount is not None:
        ecount = ae[0]
        del ae[0]
    else:
        ecount = ec
    cg.ep.count = ecount

    if rec is not None:
        cg.ep.rec = ae[0]
        del ae[0]

    if drec is not None:
        cg.ep.drec = ae[0]
        del ae[0]

    cg = Graph(cg, vorder=br)

    cg.add_vertex(B - cg.num_vertices())
    return cg

def get_entropy_args(kargs):
    kargs = kargs.copy()
    args = DictState(kargs)
    deg_dl_kind = args.degree_dl_kind
    del kargs["degree_dl_kind"]
    if deg_dl_kind == "entropy":
        kind = libinference.deg_dl_kind.ent
    elif deg_dl_kind == "uniform":
        kind = libinference.deg_dl_kind.uniform
    elif deg_dl_kind == "distributed":
        kind = libinference.deg_dl_kind.dist
    ea = libinference.entropy_args()
    ea.exact = args.exact
    ea.dense = args.dense
    ea.multigraph = args.multigraph
    ea.adjacency = args.adjacency
    del kargs["exact"]
    del kargs["dense"]
    del kargs["multigraph"]
    del kargs["adjacency"]
    if args.dl:
        ea.partition_dl = args.partition_dl
        ea.degree_dl = args.degree_dl
        ea.edges_dl = args.edges_dl
    else:
        ea.partition_dl = False
        ea.degree_dl = False
        ea.edges_dl = False
    ea.degree_dl_kind = kind
    del kargs["dl"]
    del kargs["partition_dl"]
    del kargs["degree_dl"]
    del kargs["edges_dl"]
    extract_arg(kargs, "callback", None)
    if len(kargs) > 0:
        raise ValueError("unrecognized entropy arguments: " +
                         str(list(kargs.keys())))
    return ea

_q_cache_max_n = 10000
def init_q_cache(max_n=None):
    if max_n is None:
        max_n = _q_cache_max_n
    libinference.init_q_cache(min(_q_cache_max_n, max_n))

class BlockState(object):
    r"""The stochastic block model state of a given graph.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        Graph to be modelled.
    eweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Edge multiplicities (for multigraphs or block graphs).
    vweight : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex multiplicities (for block graphs).
    rec : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Real-valued edge covariates.
    rec_type : `"positive"`, `"signed"` or `None` (optional, default: ``None``)
        Type of edge covariates. If not specified, it will be guessed from
        ``rec``.
    rec_params : ``dict`` (optional, default: ``{}``)
        Model hyperparameters for real-valued covariates. This should be a
        ``dict`` with keys in the list ``["alpha", "beta"]`` if ``rec_type ==
        positive`` or ``["m0", "k0", "v0". "nu0"]`` if ``rec_type == signed``.
    b : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Initial block labels on the vertices. If not supplied, it will be
        randomly sampled.
    B : ``int`` (optional, default: ``None``)
        Number of blocks (or vertex groups). If not supplied it will be obtained
        from the parameter ``b``.
    clabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Constraint labels on the vertices. If supplied, vertices with different
        label values will not be clustered in the same group.
    pclabel : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Partition constraint labels on the vertices. This has the same
        interpretation as ``clabel``, but will be used to compute the partition
        description length.
    deg_corr : ``bool`` (optional, default: ``True``)
        If ``True``, the degree-corrected version of the blockmodel ensemble will
        be assumed, otherwise the traditional variant will be used.
    allow_empty : ``bool`` (optional, default: ``True``)
        If ``True``, partition description length computed will allow for empty
        groups.
    max_BE : ``int`` (optional, default: ``1000``)
        If the number of blocks exceeds this value, a sparse matrix is used for
        the block graph. Otherwise a dense matrix will be used.
    """

    def __init__(self, g, eweight=None, vweight=None, rec=None, rec_type=None,
                 rec_params={}, b=None, B=None, clabel=None, pclabel=None,
                 deg_corr=True, allow_empty=False, max_BE=1000, **kwargs):
        kwargs = kwargs.copy()

        # initialize weights to unity, if necessary
        if eweight == "unity":
            eweight = g.new_ep("int", 1)
        elif eweight is None or isinstance(eweight, libinference.unity_eprop_t):
            eweight = libinference.unity_eprop_t()
        elif eweight.value_type() != "int32_t":
            eweight = g.own_property(eweight.copy(value_type="int32_t"))
        else:
            eweight = g.own_property(eweight)
        if vweight == "unity":
            vweight = g.new_vp("int", 1)
        elif vweight is None or isinstance(vweight, libinference.unity_vprop_t):
            vweight = libinference.unity_vprop_t()
        elif vweight.value_type() != "int32_t":
            vweight = g.own_property(vweight.copy(value_type="int32_t"))
        else:
            vweight = g.own_property(vweight)
        self.eweight = eweight
        self.vweight = vweight

        is_edge_weighted = not isinstance(eweight, libinference.unity_eprop_t)
        is_vertex_weighted = not isinstance(vweight, libinference.unity_vprop_t)
        if not is_edge_weighted and is_vertex_weighted:
            self.eweight = g.new_ep("int", 1)
        if not is_vertex_weighted and is_edge_weighted:
            self.vweight = g.new_vp("int", 1)
        self.is_weighted = is_edge_weighted or is_vertex_weighted

        # configure the main graph and block model parameters
        self.g = g

        self.deg_corr = deg_corr
        self.overlap = False
        self.degs = extract_arg(kwargs, "degs", libinference.simple_degs_t())
        if self.degs is None:
            self.degs = libinference.simple_degs_t()

        # ensure we have at most as many blocks as nodes
        if B is not None and b is None:
            B = min(B, self.g.num_vertices())

        if b is None:
            # create a random partition into B blocks.
            if B is None:
                raise ValueError("either 'b' or 'B' must be specified")
            B = min(B, self.g.num_vertices())
            ba = random.randint(0, B, self.g.num_vertices())
            ba[:B] = arange(B)        # avoid empty blocks
            if B < self.g.num_vertices():
                random.shuffle(ba)
            b = g.new_vp("int")
            b.fa = ba
            self.b = b
        else:
            # if a partition is available, we will incorporate it.
            if isinstance(b, numpy.ndarray):
                self.b = g.new_vp("int")
                self.b.fa = b
            else:
                self.b = b = g.own_property(b.copy(value_type="int32_t"))
            if B is None:
                B = int(self.b.fa.max()) + 1

        if self.b.fa.max() >= B:
            raise ValueError("Maximum value of b is larger or equal to B! (%d vs %d)" %
                             (self.b.fa.max(), B))

        if rec is None:
            self.rec = self.g.new_ep("double")
        else:
            self.rec = rec.copy("double")
        self.drec = extract_arg(kwargs, "drec", None)
        if self.drec is None:
            self.drec = self.rec.copy()
            self.drec.fa **= 2
        else:
            self.drec = self.drec.copy()

        # Construct block-graph
        self.bg = get_block_graph(g, B, self.b, self.vweight, self.eweight,
                                  rec=self.rec, drec=self.drec)
        self.bg.set_fast_edge_removal()

        self.mrs = self.bg.ep["count"]
        self.wr = self.bg.vp["count"]

        self.mrp = self.bg.degree_property_map("out", weight=self.mrs)

        if g.is_directed():
            self.mrm = self.bg.degree_property_map("in", weight=self.mrs)
        else:
            self.mrm = self.mrp

        self.B = B

        if pclabel is not None:
            if isinstance(pclabel, PropertyMap):
                self.pclabel = self.g.own_property(pclabel).copy("int")
            else:
                self.pclabel = self.g.new_vp("int")
                self.pclabel.fa = pclabel
        else:
            self.pclabel = self.g.new_vp("int")

        if clabel is not None:
            if isinstance(clabel, PropertyMap):
                self.clabel = self.g.own_property(clabel).copy("int")
            else:
                self.clabel = self.g.new_vp("int")
                self.clabel.fa = clabel
        elif self.pclabel.fa.max() > 0:
            self.clabel = self.pclabel
        else:
            self.clabel = self.g.new_vp("int")

        if not self._check_clabel():
            raise ValueError("provided clabel is inconsistent with node partition")

        self.bclabel = self.get_bclabel()

        self.max_BE = max_BE

        self.use_hash = self.B > self.max_BE

        self.ignore_degrees = extract_arg(kwargs, "ignore_degrees", None)
        if self.ignore_degrees is None:
            self.ignore_degrees = g.new_vp("bool", False)
        else:
            self.ignore_degrees = self.ignore_degrees.copy("bool")

        self.merge_map = extract_arg(kwargs, "merge_map",
                                     self.g.vertex_index.copy("int"))

        self.candidate_blocks = Vector_size_t()
        self.candidate_pos = self.bg.new_vp("int")
        self.empty_blocks = Vector_size_t()
        self.empty_pos = self.bg.new_vp("int")

        if rec is None:
            self.rec_type = libinference.rec_type.none
        elif rec_type is None:
            if self.rec.fa.min() > 0:
                self.rec_type = libinference.rec_type.positive
            else:
                self.rec_type = libinference.rec_type.signed
        elif rec_type == "positive":
            self.rec_type = libinference.rec_type.positive
        elif rec_type == "signed":
            self.rec_type = libinference.rec_type.signed
        elif rec_type == "delta_t":
            self.rec_type = libinference.rec_type.delta_t
        else:
            self.rec_type = rec_type

        if self.rec_type == libinference.rec_type.delta_t: # waiting times
            self.brecsum = self.bg.degree_property_map("out", self.bg.ep.rec)
            tokens = self.ignore_degrees.fa == 0
            token_groups = bincount(self.b.fa[tokens]) > 0
            token_groups.resize(self.bg.num_vertices())
            self.brecsum.a[token_groups] = 0
        else:
            self.brecsum = self.bg.new_vp("double")

        self.brec = self.bg.ep.rec
        self.bdrec = self.bg.ep.drec

        self.rec_params = dict(m0=self.rec.fa.mean(), k0=1,
                               v0=self.rec.fa.std() ** 2, nu0=3)
        if self.is_weighted:
            idx = self.eweight.fa > 0
            self.rec_params.update(dict(alpha=1, beta=self.rec.fa[idx].mean()))
        else:
            self.rec_params.update(dict(alpha=1, beta=self.rec.fa.mean()))
        self.rec_params.update(rec_params)
        self.__dict__.update(self.rec_params)

        self.allow_empty = allow_empty
        self._abg = self.bg._get_any()
        self._avweight = self.vweight._get_any()
        self._aeweight = self.eweight._get_any()
        self._state = libinference.make_block_state(self, _get_rng())

        if deg_corr:
            init_q_cache(max(self.get_E(), self.get_N()) + 1)

        self._entropy_args = dict(adjacency=True, dl=True, partition_dl=True,
                                  degree_dl=True, degree_dl_kind="distributed",
                                  edges_dl=True, dense=False, multigraph=True,
                                  exact=True)

        if len(kwargs) > 0:
            raise ValueError("unrecognized keyword arguments: " +
                             str(list(kwargs.keys())))


    def __repr__(self):
        return "<BlockState object with %d blocks (%d nonempty),%s%s for graph %s, at 0x%x>" % \
            (self.B, self.get_nonempty_B(),
             " degree-corrected," if self.deg_corr else "",
             ((" with %s real-typed edge covariates," %
               ("positive" if self.rec_type == libinference.rec_type.positive
                else "signed"))
              if self.rec_type != libinference.rec_type.none else ""),
             str(self.g), id(self))

    def __copy__(self):
        return self.copy()

    def __deepcopy__(self, memo):
        g = copy.deepcopy(self.g, memo)
        return self.copy(g=g)

    def copy(self, g=None, eweight=None, vweight=None, b=None, B=None,
             deg_corr=None, clabel=None, overlap=False, pclabel=None,
             max_BE=None, **kwargs):
        r"""Copies the block state. The parameters override the state properties, and
         have the same meaning as in the constructor."""

        if not overlap:
            state = BlockState(self.g if g is None else g,
                               eweight=self.eweight if eweight is None else eweight,
                               vweight=self.vweight if vweight is None else vweight,
                               b=self.b.copy() if b is None else b,
                               B=(self.B if b is None else None) if B is None else B,
                               clabel=self.clabel if clabel is None else clabel,
                               pclabel=self.pclabel if pclabel is None else pclabel,
                               deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                               max_BE=self.max_BE if max_BE is None else max_BE,
                               ignore_degrees=kwargs.pop("ignore_degrees",
                                                         self.ignore_degrees),
                               degs=self.degs.copy(),
                               merge_map=kwargs.get("merge_map",
                                                    self.merge_map.copy()),
                               rec=kwargs.get("rec", self.rec),
                               drec=kwargs.get("drec", self.drec),
                               rec_type=kwargs.get("rec_type", self.rec_type),
                               rec_params=kwargs.get("rec_params",
                                                     self.rec_params),
                               allow_empty=kwargs.get("allow_empty",
                                                      self.allow_empty),
                               **dmask(kwargs,
                                       ["ignore_degrees", "merge_map", "rec",
                                        "rec_type", "rec_params", "drec",
                                        "allow_empty"]))
        else:
            state = OverlapBlockState(self.g if g is None else g,
                                      b=self.b.copy() if b is None else b,
                                      B=(self.B if b is None else None) if B is None else B,
                                      rec=kwargs.get("rec", self.rec),
                                      drec=kwargs.get("drec", self.drec),
                                      rec_type=kwargs.get("rec_type",
                                                          self.rec_type),
                                      rec_params=kwargs.get("rec_params",
                                                            self.rec_params),
                                      clabel=self.clabel if clabel is None else clabel,
                                      pclabel=self.pclabel if pclabel is None else pclabel,
                                      deg_corr=self.deg_corr if deg_corr is None else deg_corr,
                                      allow_empty=kwargs.get("allow_empty",
                                                             self.allow_empty),
                                      max_BE=self.max_BE if max_BE is None else max_BE,
                                      **dmask(kwargs, ["rec", "rec_type",
                                                       "rec_params", "drec",
                                                       "allow_empty"]))
        return state


    def __getstate__(self):
        state = dict(g=self.g,
                     eweight=self.eweight if self.is_weighted else None,
                     vweight=self.vweight if self.is_weighted else None,
                     b=self.b,
                     B=self.B,
                     clabel=self.clabel,
                     pclabel=self.pclabel,
                     deg_corr=self.deg_corr,
                     allow_empty=self.allow_empty,
                     max_BE=self.max_BE,
                     rec=self.rec if self.rec_type != libinference.rec_type.none else None,
                     drec=self.drec if self.rec_type == libinference.rec_type.signed else None,
                     rec_type=int(self.rec_type),
                     rec_params=self.rec_params,
                     ignore_degrees=self.ignore_degrees.copy("int"),
                     degs=self.degs if not isinstance(self.degs,
                                                      libinference.simple_degs_t) else None,
                     merge_map=self.merge_map)
        return state

    def __setstate__(self, state):
        conv_pickle_state(state)
        self.__init__(**state)

    def get_block_state(self, b=None, vweight=False, **kwargs):
        r"""Returns a :class:`~graph_tool.community.BlockState`` corresponding to the
        block graph (i.e. the blocks of the current state become the nodes). The
        parameters have the same meaning as the in the constructor. If ``vweight
        == True`` the nodes of the block state are weighted with the node
        counts."""

        deg_corr = kwargs.get("deg_corr", self.deg_corr)
        copy_bg = extract_arg(kwargs, "copy_bg", True)

        if deg_corr and vweight:
            if isinstance(self.degs, libinference.simple_degs_t):
                degs = libinference.get_block_degs(self.g._Graph__graph,
                                                   _prop("v", self.g, self.b),
                                                   self.eweight._get_any())
            else:
                degs = libinference.get_weighted_block_degs(self.g._Graph__graph,
                                                            self.degs,
                                                            _prop("v", self.g,
                                                                  self.b))
        else:
            degs = libinference.simple_degs_t()

        if copy_bg:
            bg = self.bg.copy()
            eweight = bg.own_property(self.mrs.copy())
        else:
            bg = self.bg
            eweight = self.mrs
        if vweight == "nonempty":
            vweight = bg.new_vp("int", self.wr.a > 0)
        elif vweight == "unity":
            vweight = bg.new_vp("int", 1)
        elif vweight == True:
            if copy_bg:
                vweight = bg.own_property(self.wr.copy())
            else:
                vweight = self.wr
        else:
            vweight = None
        state = BlockState(bg,
                           eweight=eweight,
                           vweight=vweight,
                           b=bg.vertex_index.copy("int") if b is None else b,
                           deg_corr=deg_corr,
                           allow_empty=kwargs.get("allow_empty",
                                                  self.allow_empty),
                           degs=degs,
                           rec_type=kwargs.get("rec_type", self.rec_type if vweight else None),
                           rec=kwargs.get("rec", bg.ep.rec if vweight else None),
                           drec=kwargs.get("drec", bg.ep.drec if vweight else None),
                           rec_params=kwargs.get("rec_params", self.rec_params),
                           max_BE=self.max_BE,
                           **dmask(kwargs, ["allow_empty", "rec_type",
                                            "rec_params", "rec", "drec",
                                            "deg_corr"]))
        return state

    def get_E(self):
        r"Returns the total number of edges."
        return int(self.eweight.fa.sum()) if self.is_weighted else self.g.num_edges()

    def get_N(self):
        r"Returns the total number of nodes."
        return int(self.vweight.fa.sum()) if self.is_weighted else self.g.num_vertices()

    def get_B(self):
        r"Returns the total number of blocks."
        return self.bg.num_vertices()

    def get_nonempty_B(self):
        r"Returns the total number of nonempty blocks."
        return int((self.wr.a > 0).sum())

    def get_bclabel(self):
        r"""Returns a :class:`~graph_tool.PropertyMap` corresponding to constraint
        labels for the block graph."""

        bclabel = self.bg.new_vertex_property("int")
        reverse_map(self.b, bclabel)
        pmap(bclabel, self.clabel)
        return bclabel

    def get_bpclabel(self):
        r"""Returns a :class:`~graph_tool.PropertyMap`` corresponding to partition
        constraint labels for the block graph."""

        bclabel = self.bg.new_vertex_property("int")
        reverse_map(self.b, bclabel)
        pmap(bclabel, self.pclabel)
        return bclabel

    def _check_clabel(self):
        b = self.b.fa + self.clabel.fa * self.B
        b2 = self.b.fa.copy()
        continuous_map(b)
        continuous_map(b2)
        return (b == b2).all()

    def _couple_state(self, state, entropy_args):
        if state is None:
            self._state.decouple_state()
        else:
            self._state.couple_state(state._state, entropy_args)

    def get_blocks(self):
        r"""Returns the property map which contains the block labels for each vertex."""
        return self.b

    def set_blocks(self, b):
        r"""Sets the internal partition of the state."""
        if b.value_type() != "int32_t":
            b = b.copy("int32_t")
        self._state.set_partition(_prop("v", self.g, b))

    def get_bg(self):
        r"""Returns the block graph."""
        return self.bg

    def get_ers(self):
        r"""Returns the edge property map of the block graph which contains the
        :math:`e_{rs}` matrix entries.  For undirected graphs, the diagonal
        values (self-loops) contain :math:`e_{rr}/2`."""
        return self.mrs

    def get_er(self):
        r"""Returns the vertex property map of the block graph which contains the number
        :math:`e_r` of half-edges incident on block :math:`r`. If the graph is
        directed, a pair of property maps is returned, with the number of
        out-edges :math:`e^+_r` and in-edges :math:`e^-_r`, respectively."""
        if self.bg.is_directed():
            return self.mrp, self.mrm
        else:
            return self.mrp

    def get_nr(self):
        r"""Returns the vertex property map of the block graph which contains the block
        sizes :math:`n_r`."""
        return self.wr

    def entropy(self, adjacency=True, dl=True, partition_dl=True,
                degree_dl=True, degree_dl_kind="distributed", edges_dl=True,
                dense=False, multigraph=True, deg_entropy=True, exact=True,
                **kwargs):
        r"""Calculate the entropy associated with the current block partition.

        Parameters
        ----------
        adjacency : ``bool`` (optional, default: ``True``)
            If ``True``, the adjacency term of the description length will be
            included.
        dl : ``bool`` (optional, default: ``True``)
            If ``True``, the description length for the parameters will be
            included.
        partition_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the partition description length
            will be included.
        degree_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the degree sequence description
            length will be included (for degree-corrected models).
        degree_dl_kind : ``str`` (optional, default: ``"distributed"``)
            This specifies the prior used for the degree sequence. It must be
            one of: ``"uniform"``, ``"distributed"`` (default) or ``"entropy"``.
        edges_dl : ``bool`` (optional, default: ``True``)
            If ``True``, and ``dl == True`` the edge matrix description length
            will be included.
        dense : ``bool`` (optional, default: ``False``)
            If ``True``, the "dense" variant of the entropy will be computed.
        multigraph : ``bool`` (optional, default: ``True``)
            If ``True``, the multigraph entropy will be used.
        deg_entropy : ``bool`` (optional, default: ``True``)
            If ``True``, the degree entropy term that is independent of the
            network partition will be included (for degree-corrected models).
        exact : ``bool`` (optional, default: ``True``)
            If ``True``, the exact expressions will be used. Otherwise,
            Stirling's factorial approximation will be used for some terms.

        Notes
        -----

        The "entropy" of the state is minus the log-likelihood of the
        microcanonical SBM, that includes the generated graph
        :math:`\boldsymbol{A}` and the model parameters :math:`\boldsymbol{\theta}`,

        .. math::

           \mathcal{S} &= - \ln P(\boldsymbol{A},\boldsymbol{\theta}) \\
                       &= - \ln P(\boldsymbol{A}|\boldsymbol{\theta}) - \ln P(\boldsymbol{\theta}).

        This value is also called the `description length
        <https://en.wikipedia.org/wiki/Minimum_description_length>`_ of the data,
        and it corresponds to the amount of information required to describe it
        (in `nats <https://en.wikipedia.org/wiki/Nat_(unit)>`_).

        For the traditional blockmodel (``deg_corr == False``), the model
        parameters are :math:`\boldsymbol{\theta} = \{\boldsymbol{e},
        \boldsymbol{b}\}`, where :math:`\boldsymbol{e}` is the matrix of edge
        counts between blocks, and :math:`\boldsymbol{b}` is the partition of the
        nodes into blocks. For the degree-corrected blockmodel (``deg_corr ==
        True``), we have an additional set of parameters, namely the degree
        sequence :math:`\boldsymbol{k}`.

        For the traditional blockmodel, the model likelihood is

        .. math::

            P(\boldsymbol{A}|\boldsymbol{e},\boldsymbol{b}) &= \frac{\prod_{r<s}e_{rs}!\prod_re_{rr}!!}{\prod_rn_r^{e_r}}\times \frac{1}{\prod_{i<j}A_{ij}!\prod_iA_{ii}!!},\\
            P(\boldsymbol{A}|\boldsymbol{e},\boldsymbol{b}) &= \frac{\prod_{rs}e_{rs}!}{\prod_rn_r^{e_r}}\times \frac{1}{\prod_{ij}A_{ij}!},

        for undirected and directed graphs, respectively, where :math:`e_{rs}`
        is the number of edges from block :math:`r` to :math:`s` (or the number
        of half-edges for the undirected case when :math:`r=s`), and :math:`n_r`
        is the number of vertices in block :math:`r` .

        For the degree-corrected variant the equivalent expressions are

        .. math::

            P(\boldsymbol{A}|\boldsymbol{e},\boldsymbol{b},\boldsymbol{k}) &= \frac{\prod_{r<s}e_{rs}!\prod_re_{rr}!!}{\prod_re_r!}\times \frac{\prod_ik_i!}{\prod_{i<j}A_{ij}!\prod_iA_{ii}!!},\\
            P(\boldsymbol{A}|\boldsymbol{e},\boldsymbol{b},\boldsymbol{k}) &= \frac{\prod_{rs}e_{rs}!}{\prod_re_r^+!\prod_re_r^-!}\times \frac{\prod_ik_i^+!\prod_ik_i^-!}{\prod_{ij}A_{ij}!},

        where :math:`e_r = \sum_se_{rs}` is the number of half-edges incident on
        block :math:`r`, and :math:`e^+_r = \sum_se_{rs}` and :math:`e^-_r =
        \sum_se_{sr}` are the numbers of out- and in-edges adjacent to block
        :math:`r`, respectively.

        If ``exact == False``, `Stirling's approximation
        <https://en.wikipedia.org/wiki/Stirling%27s_approximation>`_ is used in
        the above expression.

        If ``dense == True``, the likelihood for the non-degree-corrected model
        becomes instead

        .. math::

            P(\boldsymbol{A}|\boldsymbol{e},\boldsymbol{b})^{-1} &= \prod_{r<s}{n_rn_s\choose e_{rs}}\prod_r{{n_r\choose 2}\choose e_{rr}/2},\\
            P(\boldsymbol{A}|\boldsymbol{e},\boldsymbol{b})^{-1} &= \prod_{rs}{n_rn_s\choose e_{rs}}

        if ``multigraph == False``, otherwise we replace :math:`{n\choose
        m}\to\left(\!\!{n\choose m}\!\!\right)` above, where
        :math:`\left(\!\!{n\choose m}\!\!\right) = {n+m-1\choose m}`.  A dense
        entropy for the degree-corrected model is not available, and if
        requested will raise a :exc:`NotImplementedError`.

        If ``dl == True``, the description length :math:`\mathcal{L} = -\ln
        P(\boldsymbol{\theta})` of the model will be returned as well. The terms
        :math:`P(\boldsymbol{e})` and :math:`P(\boldsymbol{b})` are described in
        described in :func:`model_entropy`.

        For the degree-corrected model we need to specify the prior
        :math:`P(\boldsymbol{k})` for the degree sequence as well. Here there
        are three options:

        1. ``degree_dl_kind == "uniform"``

            .. math::

                P(\boldsymbol{k}|\boldsymbol{e},\boldsymbol{b}) = \prod_r\left(\!\!{n_r\choose e_r}\!\!\right)^{-1}.

        2. ``degree_dl_kind == "distributed"``

            .. math::

                P(\boldsymbol{k}|\boldsymbol{e},\boldsymbol{b}) = \prod_r\frac{\prod_k\eta_k^r!}{n_r!} \prod_r q(e_r, n_r)^{-1}

            with :math:`\eta_k^r` being the number of nodes with degree
            :math:`k` in group :math:`r`, and :math:`q(n,m)` being the number of
            `partitions
            <https://en.wikipedia.org/wiki/Partition_(number_theory)>`_ of
            integer :math:`n` into at most :math:`m` parts.

        3. ``degree_dl_kind == "entropy"``

            .. math::

                P(\boldsymbol{k}|\boldsymbol{e},\boldsymbol{b}) \approx \prod_r\exp\left(-n_rH(\boldsymbol{k}_r)\right)

            where :math:`H(\boldsymbol{k}_r) = -\sum_kp_r(k)\ln p_r(k)` is the
            entropy of the degree distribution inside block :math:`r`.

            Note that, differently from the other two choices, this represents
            only an approximation of the description length. It is meant to be
            used only for comparison purposes, and should be avoided in practice.


        For the directed case, the above expressions are duplicated for the in-
        and out-degrees.

        """

        if _bm_test() and kwargs.get("test", True):
            args = dict(**locals())
            args.update(**kwargs)
            del args["self"]
            del args["kwargs"]

        kwargs = kwargs.copy()

        E = self.get_E()
        N = self.get_N()

        S = 0
        if adjacency:
            S += self._state.entropy(dense, multigraph, deg_entropy, exact)

            if not dense and not exact:
                if multigraph:
                    S -= E
                else:
                    S += E

        if dl:
            if partition_dl:
                S += self._state.get_partition_dl()

            if edges_dl:
                if self.allow_empty:
                    actual_B = self.B
                else:
                    actual_B = (self.wr.a > 0).sum()
                S += model_entropy(actual_B, N, E,
                                   directed=self.g.is_directed(), nr=False)

            if self.deg_corr and degree_dl:
                if degree_dl_kind == "entropy":
                    kind = libinference.deg_dl_kind.ent
                elif degree_dl_kind == "uniform":
                    kind = libinference.deg_dl_kind.uniform
                elif degree_dl_kind == "distributed":
                    kind = libinference.deg_dl_kind.dist
                S += self._state.get_deg_dl(kind)

        callback = extract_arg(kwargs, "callback", None)
        if callback is not None:
            S += callback(self)

        if extract_arg(kwargs, "test", True) and _bm_test():
            assert not isnan(S) and not isinf(S), \
                "invalid entropy %g (%s) " % (S, str(args))

            args["test"] = False
            state_copy = self.copy()
            Salt = state_copy.entropy(**args)
            assert abs(S - Salt) < 1e-6, \
                "entropy discrepancy after copying (%g %g)" % (S, Salt)

        if len(kwargs) > 0:
            raise ValueError("unrecognized keyword arguments: " +
                             str(list(kwargs.keys())))
        return S

    def get_matrix(self):
        r"""Returns the block matrix (as a sparse :class:`~scipy.sparse.csr_matrix`),
        which contains the number of edges between each block pair.

        .. warning::

           This corresponds to the adjacency matrix of the block graph, which by
           convention includes twice the amount of edges in the diagonal entries
           if the graph is undirected.

        Examples
        --------

        .. testsetup:: get_matrix

           gt.seed_rng(42)
           np.random.seed(42)
           from pylab import *

        .. doctest:: get_matrix

           >>> g = gt.collection.data["polbooks"]
           >>> state = gt.BlockState(g, B=5, deg_corr=True)
           >>> state.mcmc_sweep(niter=1000)
           (...)
           >>> m = state.get_matrix()
           >>> figure()
           <...>
           >>> matshow(m.todense())
           <...>
           >>> savefig("bloc_mat.pdf")

        .. testcleanup:: get_matrix

           savefig("bloc_mat.png")

        .. figure:: bloc_mat.*
           :align: center

           A  5x5 block matrix.

       """

        return adjacency(self.bg, weight=self.mrs)

    def virtual_vertex_move(self, v, s, **kwargs):
        r"""Computes the entropy difference if vertex ``v`` is moved to block ``s``. The
        remaining parameters are the same as in
        :meth:`graph_tool.BlockState.entropy`."""
        return self._state.virtual_move(int(v), self.b[v], s,
                                        get_entropy_args(overlay(self._entropy_args,
                                                                 **kwargs)))

    def move_vertex(self, v, s):
        r"""Move vertex ``v`` to block ``s``.

        This optionally accepts a list of vertices and blocks to move
        simultaneously.
        """
        if not isinstance(v, collections.Iterable):
            self._state.move_vertex(int(v), s)
        else:
            self._state.move_vertices(numpy.asarray(v, dtype="uint64"),
                                      numpy.asarray(s, dtype="uint64"))

    def remove_vertex(self, v):
        r"""Remove vertex ``v`` from its current group.

        This optionally accepts a list of vertices to remove.

        .. warning::

           This will leave the state in an inconsistent state before the vertex
           is returned to some other group, or if the same vertex is removed
           twice.
        """
        if isinstance(v, collections.Iterable):
            if not isinstance(v, numpy.ndarray):
                v = list(v)
            self._state.remove_vertices(numpy.asarray(v, dtype="uint64"))
        else:
            self._state.remove_vertex(int(v))

    def add_vertex(self, v, r):
        r"""Add vertex ``v`` to block ``r``.

        This optionally accepts a list of vertices and blocks to add.

        .. warning::

           This can leave the state in an inconsistent state if a vertex is
           added twice to the same group.
        """
        if isinstance(v, collections.Iterable):
            if not isinstance(v, numpy.ndarray):
                v = list(v)
            if not isinstance(r, numpy.ndarray):
                r = list(r)
            self._state.add_vertices(numpy.asarray(v, dtype="uint64"),
                                     numpy.asarray(r, dtype="uint64"))
        else:
            self._state.add_vertex(int(v), r)

    def merge_vertices(self, u, v):
        r"""Merge vertex ``u`` into ``v``.

        .. warning::

           This modifies the underlying graph.
        """
        self.move_vertex(u, self.b[v])
        self._state.merge_vertices(int(u), int(v))

    def sample_vertex_move(self, v, c=1.):
        r"""Sample block membership proposal of vertex ``v`` according to real-valued
        sampling parameter ``c``: For :math:`c\to 0` the blocks are sampled
        according to the local neighbourhood and their connections; for
        :math:`c\to\infty` the blocks are sampled randomly.
        """
        return self._state.sample_block(int(v), c, _get_rng())

    def get_move_prob(self, v, s, c=1., reverse=False):
        r"""Compute the probability of a move proposal for vertex ``v`` to block ``s``
        according to sampling parameter ``c``, as obtained with
        :meth:`graph_tool.inference.BlockState.sample_vertex_move`. If
        ``reverse == True``, the reverse probability of moving the node back
        from block ``s`` to its current one is obtained.
        """
        if not reverse:
            return self._state.get_move_prob(int(v), self.b[v], s, c, False)
        else:
            return self._state.get_move_prob(int(v), s, self.b[v], c, True)

    def get_edges_prob(self, edge_list, missing=True, entropy_args={}):
        """Compute the unnormalized log-probability of the missing (or spurious if
        ``missing=False``) edges given by ``edge_list`` (a list of ``(source,
        target)`` tuples, or :meth:`~graph_tool.Edge` instances). The values in
        ``entropy_args`` are passed to :meth:`graph_tool.BlockState.entropy()`
        to calculate the log-probability.

        """
        pos = {}
        for u, v in edge_list:
            pos[u] = self.b[u]
            pos[v] = self.b[v]

        Si = self.entropy(**entropy_args)

        self.remove_vertex(pos.keys())

        try:
            if missing:
                new_es = []
                for u, v in edge_list:
                    if not self.is_weighted:
                        e = self.g.add_edge(u, v)
                    else:
                        e = self.g.edge(u, v)
                        if e is None:
                            e = self.g.add_edge(u, v)
                            self.eweight[e] = 0
                        self.eweight[e] += 1
                    new_es.append(e)
            else:
                old_es = []
                for e in edge_list:
                    u, v = e
                    if isinstance(e, tuple):
                        e = self.g.edge(u, v)
                        if e is None:
                            raise ValueError("edge not found: (%d, %d)" % (int(u),
                                                                           int(v)))

                    if self.is_weighted:
                        self.eweight[e] -= 1
                        if self.eweight[e] == 0:
                            self.g.remove_edge(e)
                    else:
                        self.g.remove_edge(e)
                    old_es.append((u, v))

            self.add_vertex(pos.keys(), pos.values())

            Sf = self.entropy(**entropy_args)

            self.remove_vertex(pos.keys())

        finally:
            if missing:
                if self.is_weighted:
                    for e in reversed(new_es):
                        self.eweight[e] -= 1
                        if self.eweight[e] == 0:
                            self.g.remove_edge(e)
                else:
                    for e in reversed(new_es):
                        self.g.remove_edge(e)
            else:
                for u, v in old_es:
                    if self.is_weighted:
                        e = self.g.edge(u, v)
                        if e is None:
                            e = self.g.add_edge(u, v)
                            self.eweight[e] = 0
                        self.eweight[e] += 1
                    else:
                        self.g.add_edge(u, v)

            self.add_vertex(pos.keys(), pos.values())

        L = Si - Sf

        if _bm_test():
            state = self.copy()
            set_test(False)
            L_alt = state.get_edges_prob(edge_list, missing=missing,
                                         entropy_args=entropy_args)
            set_test(True)
            assert abs(L - L_alt) < 1e-6, \
                "inconsistent missing=%s edge probability (%g, %g): %s, %s" % \
                (str(missing), L, L_alt,  str(entropy_args), str(edge_list))

        return L

    def _mcmc_sweep_dispatch(self, mcmc_state):
        return libinference.mcmc_sweep(mcmc_state, self._state,
                                       _get_rng())

    def mcmc_sweep(self, beta=1., c=1., niter=1, entropy_args={},
                   allow_vacate=True, sequential=True, parallel=False,
                   vertices=None, verbose=False, **kwargs):
        r"""Perform ``niter`` sweeps of a Metropolis-Hastings acceptance-rejection
        sampling MCMC to sample network partitions.

        Parameters
        ----------
        beta : ``float`` (optional, default: ``1.``)
            Inverse temperature.
        c : ``float`` (optional, default: ``1.``)
            Sampling parameter ``c`` for move proposals: For :math:`c\to 0` the
            blocks are sampled according to the local neighbourhood of a given
            node and their block connections; for :math:`c\to\infty` the blocks
            are sampled randomly. Note that only for :math:`c > 0` the MCMC is
            guaranteed to be ergodic.
        niter : ``int`` (optional, default: ``1``)
            Number of sweeps to perform. During each sweep, a move attempt is
            made for each node.
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        allow_vacate : ``bool`` (optional, default: ``True``)
            Allow groups to be vacated.
        sequential : ``bool`` (optional, default: ``True``)
            If ``sequential == True`` each vertex move attempt is made
            sequentially, where vertices are visited in random order. Otherwise
            the moves are attempted by sampling vertices randomly, so that the
            same vertex can be moved more than once, before other vertices had
            the chance to move.
        parallel : ``bool`` (optional, default: ``False``)
            If ``parallel == True``, vertex movements are attempted in parallel.

            .. warning::

               If ``parallel == True``, the asymptotic exactness of the MCMC
               sampling is not guaranteed.
        vertices : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of vertices which will be
            moved. Otherwise, all vertices will.
        verbose : ``bool`` (optional, default: ``False``)
            If ``verbose == True``, detailed information will be displayed.

        Returns
        -------
        dS : ``float``
            Entropy difference after the sweeps.
        nmoves : ``int``
            Number of vertices moved.

        Notes
        -----
        This algorithm has an :math:`O(E)` complexity, where :math:`E` is the
        number of edges (independent of the number of blocks).

        References
        ----------
        .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and
           greedy heuristic for the inference of stochastic block models", Phys.
           Rev. E 89, 012804 (2014), :doi:`10.1103/PhysRevE.89.012804`,
           :arxiv:`1310.4378`
        """

        mcmc_state = DictState(locals())
        entropy_args = overlay(self._entropy_args, **entropy_args)
        if (_bm_test() and entropy_args["multigraph"] and
            not entropy_args["dense"] and
            hasattr(self, "degs") and
            not isinstance(self.degs, libinference.simple_degs_t)):
            entropy_args["multigraph"] = False
        mcmc_state.entropy_args = get_entropy_args(entropy_args)
        mcmc_state.vlist = Vector_size_t()
        if vertices is None:
            vertices = self.g.vertex_index.copy().fa
            if self.is_weighted:
                # ignore vertices with zero weight
                vw = self.vweight.fa
                vertices = vertices[vw > 0]
        mcmc_state.vlist.resize(len(vertices))
        mcmc_state.vlist.a = vertices
        mcmc_state.E = self.get_E()
        mcmc_state.state = self._state

        disable_callback_test = extract_arg(kwargs, "disable_callback_test", False)
        if _bm_test():
            assert self._check_clabel(), "invalid clabel before sweep"
            if disable_callback_test:
                Si = self.entropy(**dmask(entropy_args, ["callback"]))
            else:
                Si = self.entropy(**entropy_args)

        dS, nmoves = self._mcmc_sweep_dispatch(mcmc_state)

        if _bm_test():
            assert self._check_clabel(), "invalid clabel after sweep"
            if disable_callback_test:
                Sf = self.entropy(**dmask(entropy_args, ["callback"]))
            else:
                Sf = self.entropy(**entropy_args)
            assert abs(dS - (Sf - Si)) < 1e-6, \
                "inconsistent entropy delta %g (%g): %s" % (dS, Sf - Si,
                                                            str(entropy_args))

        if len(kwargs) > 0:
            raise ValueError("unrecognized keyword arguments: " +
                             str(list(kwargs.keys())))

        return dS, nmoves

    def _gibbs_sweep_dispatch(self, gibbs_state):
        return libinference.gibbs_sweep(gibbs_state, self._state,
                                        _get_rng())

    def gibbs_sweep(self, beta=1., niter=1, entropy_args={}, allow_vacate=True,
                    sequential=True, parallel=False, vertices=None,
                    verbose=False, **kwargs):
        r"""Perform ``niter`` sweeps of a rejection-free Gibbs sampling MCMC
        to sample network partitions.

        Parameters
        ----------
        beta : ``float`` (optional, default: ``1.``)
            Inverse temperature.
        niter : ``int`` (optional, default: ``1``)
            Number of sweeps to perform. During each sweep, a move attempt is
            made for each node.
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        allow_vacate : ``bool`` (optional, default: ``True``)
            Allow groups to be vacated.
        sequential : ``bool`` (optional, default: ``True``)
            If ``sequential == True`` each vertex move attempt is made
            sequentially, where vertices are visited in random order. Otherwise
            the moves are attempted by sampling vertices randomly, so that the
            same vertex can be moved more than once, before other vertices had
            the chance to move.
        parallel : ``bool`` (optional, default: ``False``)
            If ``parallel == True``, vertex movements are attempted in parallel.

            .. warning::

               If ``parallel == True``, the asymptotic exactness of the MCMC
               sampling is not guaranteed.
        vertices : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of vertices which will be
            moved. Otherwise, all vertices will.
        verbose : ``bool`` (optional, default: ``False``)
            If ``verbose == True``, detailed information will be displayed.

        Returns
        -------
        dS : ``float``
            Entropy difference after the sweeps.
        nmoves : ``int``
            Number of vertices moved.

        Notes
        -----
        This algorithm has an :math:`O(E\times B)` complexity, where :math:`B`
        is the number of blocks, and :math:`E` is the number of edges.

        """

        gibbs_state = DictState(locals())
        entropy_args = overlay(self._entropy_args, **entropy_args)
        gibbs_state.entropy_args = get_entropy_args(entropy_args)
        gibbs_state.vlist = Vector_size_t()
        if vertices is None:
            vertices = self.g.vertex_index.copy().fa
            if self.is_weighted:
                # ignore vertices with zero weight
                vw = self.vweight.fa
                vertices = vertices[vw > 0]
        gibbs_state.vlist.resize(len(vertices))
        gibbs_state.vlist.a = vertices
        gibbs_state.E = self.get_E()
        gibbs_state.state = self._state

        disable_callback_test = extract_arg(kwargs, "disable_callback_test", False)
        if _bm_test():
            assert self._check_clabel(), "invalid clabel before sweep"
            if disable_callback_test:
                Si = self.entropy(**dmask(entropy_args, ["callback"]))
            else:
                Si = self.entropy(**entropy_args)

        dS, nmoves = self._gibbs_sweep_dispatch(gibbs_state)

        if _bm_test():
            assert self._check_clabel(), "invalid clabel after sweep"
            if disable_callback_test:
                Sf = self.entropy(**dmask(entropy_args, ["callback"]))
            else:
                Sf = self.entropy(**entropy_args)
            assert abs(dS - (Sf - Si)) < 1e-6, \
                "inconsistent entropy delta %g (%g): %s" % (dS, Sf - Si,
                                                            str(entropy_args))

        if len(kwargs) > 0:
            raise ValueError("unrecognized keyword arguments: " +
                             str(list(kwargs.keys())))

        return dS, nmoves

    def _multicanonical_sweep_dispatch(self, multicanonical_state):
        return libinference.multicanonical_sweep(multicanonical_state,
                                                 self._state, _get_rng())

    def multicanonical_sweep(self, m_state, c=numpy.inf, niter=1,
                             entropy_args={}, allow_vacate=True, vertices=None,
                             verbose=False, **kwargs):
        r"""Perform ``niter`` sweeps of a non-Markovian multicanonical sampling using the
        Wang-Landau algorithm.

        Parameters
        ----------
        m_state : :class:`~graph_tool.inference.MulticanonicalState`
            :class:`~graph_tool.inference.MulticanonicalState` instance
            containing the current state of the Wang-Landau run.
        c : ``float`` (optional, default: ``numpy.inf``)
            Sampling parameter ``c`` for move proposals: For :math:`c\to 0` the
            blocks are sampled according to the local neighbourhood of a given
            node and their block connections; for :math:`c\to\infty` the blocks
            are sampled randomly. Note that only for :math:`c > 0` the MCMC is
            guaranteed to be ergodic.
        niter : ``int`` (optional, default: ``1``)
            Number of sweeps to perform. During each sweep, a move attempt is
            made for each node.
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        allow_vacate : ``bool`` (optional, default: ``True``)
            Allow groups to be vacated.
        vertices : ``list`` of ints (optional, default: ``None``)
            If provided, this should be a list of vertices which will be
            moved. Otherwise, all vertices will.
        verbose : ``bool`` (optional, default: ``False``)
            If ``verbose == True``, detailed information will be displayed.

        Returns
        -------
        dS : ``float``
            Entropy difference after the sweeps.
        nmoves : ``int``
            Number of vertices moved.

        Notes
        -----
        This algorithm has an :math:`O(E)` complexity, where :math:`E` is the
        number of edges (independent of the number of blocks).

        References
        ----------
        .. [wang-efficient-2001] Fugao Wang, D. P. Landau, "An efficient, multiple
           range random walk algorithm to calculate the density of states", Phys.
           Rev. Lett. 86, 2050 (2001), :doi:`10.1103/PhysRevLett.86.2050`,
           :arxiv:`cond-mat/0011174`

        .. [belardinelli-wang-2007] R. E. Belardinelli, V. D. Pereyra,
           "Wang-Landau algorithm: A theoretical analysis of the saturation of
           the error", J. Chem. Phys. 127, 184105 (2007),
           :doi:`10.1063/1.2803061`, :arxiv:`cond-mat/0702414`
        """

        niter *= self.g.num_vertices()
        args = dmask(locals(), ["self"])
        multi_state = DictState(args)
        entropy_args = overlay(self._entropy_args, **entropy_args)
        multi_state.entropy_args = get_entropy_args(entropy_args)
        multi_state.update(entropy_args)
        multi_state.vlist = Vector_size_t()
        if vertices is None:
            vertices = self.g.vertex_index.copy().fa
            if self.is_weighted:
                # ignore vertices with zero weight
                vw = self.vweight.fa
                vertices = vertices[vw > 0]
        multi_state.vlist.resize(len(vertices))
        multi_state.vlist.a = vertices
        multi_state.E = self.get_E()
        multi_state.S = self.entropy(**entropy_args)
        multi_state.state = self._state

        multi_state.f = m_state._f
        multi_state.time = m_state._time
        multi_state.refine = m_state._refine
        multi_state.S_min = m_state._S_min
        multi_state.S_max = m_state._S_max
        multi_state.hist = m_state._hist
        multi_state.dens = m_state._density

        S, nmoves, f, time = \
                self._multicanonical_sweep_dispatch(multi_state)

        m_state._f = f
        m_state._time = time

        disable_callback_test = extract_arg(kwargs, "disable_callback_test", False)

        if _bm_test():
            assert self._check_clabel(), "invalid clabel after sweep"
            if disable_callback_test:
                Sf = self.entropy(**dmask(entropy_args, ["callback"]))
            else:
                Sf = self.entropy(**entropy_args)
            assert abs(S - Sf) < 1e-6, \
                "inconsistent entropy after sweep %g (%g): %s" % \
                (S, Sf, str(entropy_args))

        if len(kwargs) > 0:
            raise ValueError("unrecognized keyword arguments: " +
                             str(list(kwargs.keys())))

        return S, nmoves

    def _exhaustive_sweep_dispatch(self, exhaustive_state, callback, hist):
        if callback is not None:
            return libinference.exhaustive_sweep(exhaustive_state, self._state,
                                                 callback)
        else:
            if hist is None:
                return libinference.exhaustive_sweep_iter(exhaustive_state,
                                                          self._state)
            else:
                return libinference.exhaustive_dens(exhaustive_state,
                                                    self._state, hist[0],
                                                    hist[1], hist[2])

    def exhaustive_sweep(self, entropy_args={}, callback=None, density=None,
                         vertices=None, initial_partition=None, max_iter=None):
        r"""Perform an exhaustive loop over all possible network partitions.

        Parameters
        ----------
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        callback : callable object (optional, default: ``None``)
            Function to be called for each partition, with three arguments ``(S,
            S_min, b_min)`` corresponding to the the current entropy value, the
            minimum entropy value so far, and the corresponding partition,
            respectively. If not provided, and ``hist == None`` an iterator over
            the same values will be returned instead.
        density : ``tuple`` (optional, default: ``None``)
            If provided, it should contain a tuple with values ``(S_min, S_max,
            n_bins)``, which will be used to obtain the density of states via a
            histogram of size ``n_bins``. This parameter is ignored unless
            ``callback == None``.
        vertices : iterable of ints (optional, default: ``None``)
            If provided, this should be a list of vertices which will be
            moved. Otherwise, all vertices will.
        initial_partition : iterable  of ints (optional, default: ``None``)
            If provided, this will provide the initial partition for the
            iteration.
        max_iter : ``int`` (optional, default: ``None``)
            If provided, this will limit the total number of iterations.

        Returns
        -------
        states : iterator over (S, S_min, b_min)
            If ``callback`` is ``None`` and ``hist`` is ``None``, the function
            will return an iterator over ``(S, S_min, b_min)`` corresponding to
            the the current entropy value, the minimum entropy value so far, and
            the corresponding partition, respectively.
        Ss, counts : pair of :class:`numpy.ndarray`
            If ``callback is None`` and ``hist is not None``, the function will
            return the values of each bin (``Ss``) and the state count of each
            bin (``counts``).
        b_min : :class:`~graph_tool.PropertyMap`
            If ``callback is not None`` or ``hist is not None``, the function
            will also return partition with smallest entropy.

        Notes
        -----

        This algorithm has an :math:`O(B^N)` complexity, where :math:`B` is the
        number of blocks, and :math:`N` is the number of vertices.

        """

        exhaustive_state = DictState(dict(max_iter=max_iter if max_iter is not None else 0))
        entropy_args = overlay(self._entropy_args, **entropy_args)
        exhaustive_state.entropy_args = get_entropy_args(entropy_args)
        exhaustive_state.vlist = Vector_size_t()
        if vertices is None:
            vertices = self.g.vertex_index.copy().fa
            if self.is_weighted:
                # ignore vertices with zero weight
                vw = self.vweight.fa
                vertices = vertices[vw > 0]
        if initial_partition is None:
            initial_partition = zeros(len(vertices), dtype="uint64")
        self.move_vertex(vertices, initial_partition)
        exhaustive_state.vlist.resize(len(vertices))
        exhaustive_state.vlist.a = vertices
        exhaustive_state.S = self.entropy(**entropy_args)
        exhaustive_state.state = self._state
        exhaustive_state.b_min = b_min = self.g.new_vp("int32_t")

        if density is not None:
            density = (density[0], density[1],
                       numpy.zeros(density[2], dtype="uint64"))
        if callback is not None:
            _callback = lambda S, S_min: callback(S, S_min, b_min)
        else:
            _callback = None
        ret = self._exhaustive_sweep_dispatch(exhaustive_state, _callback,
                                              density)
        if _callback is None:
            if density is None:
                return ((S, S_min, b_min) for S, S_min in ret)
            else:
                Ss = numpy.linspace(density[0], density[1], len(density[2]))
                return (Ss, density[2]), b_min
        else:
            return b_min

    def _merge_sweep_dispatch(self, merge_state):
        if not self.is_weighted:
            raise ValueError("state must be weighted to perform merges")

        if (merge_state.entropy_args.multigraph and
            not merge_state.entropy_args.dense):
            raise ValueError("can only merge multigraphs if dense == True")

        return libinference.merge_sweep(merge_state, self._state,
                                        _get_rng())

    def merge_sweep(self, nmerges=1, niter=10, entropy_args={}, parallel=True,
                    verbose=False, **kwargs):
        r"""Perform ``niter`` merge sweeps, where block nodes are progressively
        merged together in a manner that least increases the entropy.

        Parameters
        ----------
        nmerges : ``int`` (optional, default: ``1``)
            Number block nodes to merge.
        niter : ``int`` (optional, default: ``1``)
            Number of merge attempts to perform for each block node, before the
            best one is selected.
        entropy_args : ``dict`` (optional, default: ``{}``)
            Entropy arguments, with the same meaning and defaults as in
            :meth:`graph_tool.inference.BlockState.entropy`.
        parallel : ``bool`` (optional, default: ``True``)
            If ``parallel == True``, the merge candidates are obtained in
            parallel.
        verbose : ``bool`` (optional, default: ``False``)
            If ``verbose == True``, detailed information will be displayed.

        Notes
        -----

        This function should only be called for block states, obtained from
        :meth:`graph_tool.inference.BlockState.get_block_state`.

        Returns
        -------
        dS : ``float``
            Entropy difference after the sweeps.
        nmoves : ``int``
            Number of vertices merged.

        References
        ----------
        .. [peixoto-efficient-2014] Tiago P. Peixoto, "Efficient Monte Carlo and
           greedy heuristic for the inference of stochastic block models", Phys.
           Rev. E 89, 012804 (2014), :doi:`10.1103/PhysRevE.89.012804`,
           :arxiv:`1310.4378`
        """
        merge_state = DictState(locals())
        merge_state.E = self.get_E()
        merge_state.state = self._state
        entropy_args = overlay(self._entropy_args, **entropy_args)
        merge_state.entropy_args = get_entropy_args(entropy_args)

        disable_callback_test = extract_arg(kwargs, "disable_callback_test", False)
        if _bm_test():
            self._check_clabel()
            if disable_callback_test:
                Si = self.entropy(**dmask(entropy_args, ["callback"]))
            else:
                Si = self.entropy(**entropy_args)

        dS, nmoves = self._merge_sweep_dispatch(merge_state)

        if _bm_test():
            assert self._check_clabel(), "invalid clabel after sweep"
            if disable_callback_test:
                Sf = self.entropy(**dmask(entropy_args, ["callback"]))
            else:
                Sf = self.entropy(**entropy_args)
            assert abs(dS - (Sf - Si)) < 1e-6, \
                "inconsistent entropy delta %g (%g): %s" % (dS, Sf - Si,
                                                            str(entropy_args))

        return dS, nmoves

    def shrink(self, B, **kwargs):
        """Reduces the order of current state by progressively merging groups, until
        only ``B`` are left. All remaining keyword arguments are passed to
        :meth:`graph_tool.inference.BlockState.merge_sweep`.

        This function leaves the current state untouched and returns instead a
        copy with the new partition.
        """
        bstate = self.get_block_state(vweight=True,
                                      clabel=self.get_bclabel(),
                                      deg_corr=self.deg_corr)
        nB = bstate.get_nonempty_B()
        while nB > B:
            bstate.merge_sweep(nB - B, **kwargs)
            nB = bstate.get_nonempty_B()
        b = self.b.copy()
        pmap(b, bstate.merge_map)
        continuous_map(b)
        state = self.copy(b=b)
        if _bm_test():
            nB = (state.wr.a > 0).sum()
            assert nB == B, "wrong number of blocks after shrink: %d (should be %d)" % (nB, B)
        return state

    def collect_edge_marginals(self, p=None, update=1):
        r"""Collect the edge marginal histogram, which counts the number of times
        the endpoints of each node have been assigned to a given block pair.

        This should be called multiple times, e.g. after repeated runs of the
        :meth:`graph_tool.inference.BlockState.mcmc_sweep` function.

        Parameters
        ----------
        p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
            Edge property map with edge marginals to be updated.  If not
            provided, an empty histogram will be created.
        update : float (optional, default: ``1``)
            Each call increases the current count by the amount given by this
            parameter.

        Returns
        -------
        p : :class:`~graph_tool.PropertyMap`
            Edge property map with updated edge marginals.

        Examples
        --------
        .. testsetup:: collect_edge_marginals

           gt.seed_rng(42)
           np.random.seed(42)

        .. doctest:: collect_edge_marginals

           >>> g = gt.collection.data["polbooks"]
           >>> state = gt.BlockState(g, B=4, deg_corr=True)
           >>> pe = None
           >>> state.mcmc_sweep(niter=1000)   # remove part of the transient
           (...)
           >>> for i in range(1000):
           ...     ds, nmoves = state.mcmc_sweep(niter=10)
           ...     pe = state.collect_edge_marginals(pe)
           >>> gt.bethe_entropy(g, pe)[0]
           6.65429696440...
        """

        if p is None:
            p = self.g.new_ep("python::object",
                              vals=[libinference.BlockPairHist()
                                    for i in range(self.g.num_edges())])

        libinference.edge_marginals(self.g._Graph__graph,
                                    _prop("v", self.g, self.b),
                                    _prop("e", self.g, p),
                                    update)
        return p

    def collect_vertex_marginals(self, p=None, update=1):
        r"""Collect the vertex marginal histogram, which counts the number of times a
        node was assigned to a given block.

        This should be called multiple times, e.g. after repeated runs of the
        :meth:`graph_tool.inference.BlockState.mcmc_sweep` function.

        Parameters
        ----------
        p : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
            Vertex property map with vector-type values, storing the previous block
            membership counts. If not provided, an empty histogram will be created.
        update : int (optional, default: ``1``)
            Each call increases the current count by the amount given by this
            parameter.

        Returns
        -------
        p : :class:`~graph_tool.PropertyMap`
            Vertex property map with vector-type values, storing the accumulated
            block membership counts.

        Examples
        --------
        .. testsetup:: cvm

           gt.seed_rng(42)
           np.random.seed(42)

        .. doctest:: cvm

           >>> g = gt.collection.data["polbooks"]
           >>> state = gt.BlockState(g, B=4, deg_corr=True)
           >>> pv = None
           >>> state.mcmc_sweep(niter=1000)   # remove part of the transient
           (...)
           >>> for i in range(1000):
           ...     ds, nmoves = state.mcmc_sweep(niter=10)
           ...     pv = state.collect_vertex_marginals(pv)
           >>> gt.mf_entropy(g, pv)
           29.84318996...
           >>> gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie",
           ...               vertex_pie_fractions=pv, output="polbooks_blocks_soft_B4.pdf")
           <...>

        .. testcleanup:: cvm

           gt.graph_draw(g, pos=g.vp["pos"], vertex_shape="pie", vertex_pie_fractions=pv,
                         output="polbooks_blocks_soft_B4.png")

        .. figure:: polbooks_blocks_soft_B4.*
           :align: center

           "Soft" block partition of a political books network with :math:`B=4`.

        """

        if p is None:
            p = self.g.new_vp("vector<double>")

        libinference.vertex_marginals(self.g._Graph__graph,
                                      _prop("v", self.g, self.b),
                                      _prop("v", self.g, p),
                                      update)
        return p

    def collect_partition_histogram(self, h=None, update=1):
        r"""Collect a histogram of partitions.

        This should be called multiple times, e.g. after repeated runs of the
        :meth:`graph_tool.inference.BlockState.mcmc_sweep` function.

        Parameters
        ----------
        h : :class:`~graph_tool.inference.PartitionHist` (optional, default: ``None``)
            Partition histogram. If not provided, an empty histogram will be created.
        update : float (optional, default: ``1``)
            Each call increases the current count by the amount given by this
            parameter.

        Returns
        -------
        h : :class:`~graph_tool.inference.PartitionHist` (optional, default: ``None``)
            Updated Partition histogram.

        Examples
        --------
        .. testsetup:: cvm

           gt.seed_rng(42)
           np.random.seed(42)

        .. doctest:: cvm

           >>> g = gt.collection.data["polbooks"]
           >>> state = gt.BlockState(g, B=4, deg_corr=True)
           >>> ph = None
           >>> state.mcmc_sweep(niter=1000)   # remove part of the transient
           (...)
           >>> for i in range(1000):
           ...     ds, nmoves = state.mcmc_sweep(niter=10)
           ...     ph = state.collect_partition_histogram(ph)
           >>> gt.microstate_entropy(ph)
           5.215767...
        """

        if h is None:
            h = PartitionHist()
        libinference.collect_partitions(_prop("v", self.g, self.b),
                                        h, update)
        return h

    def draw(self, **kwargs):
        r"""Convenience wrapper to :func:`~graph_tool.draw.graph_draw` that
        draws the state of the graph as colors on the vertices and edges."""
        gradient = self.g.new_ep("double")
        gradient = group_vector_property([gradient])
        from graph_tool.draw import graph_draw
        return graph_draw(self.g,
                          vertex_fill_color=kwargs.get("vertex_fill_color",
                                                       self.b),
                          vertex_color=kwargs.get("vertex_color", self.b),
                          edge_gradient=kwargs.get("edge_gradient",
                                                   gradient),
                          **dmask(kwargs, ["vertex_fill_color",
                                           "vertex_color",
                                           "edge_gradient"]))

def model_entropy(B, N, E, directed=False, nr=None, allow_empty=True):
    r"""Computes the amount of information necessary for the parameters of the
    traditional blockmodel ensemble, for ``B`` blocks, ``N`` vertices, ``E``
    edges, and either a directed or undirected graph.

    This is equivalently defined as minus the log-likelihood of sampling the
    parameters from a nonparametric generative model.

    A traditional blockmodel is defined as a set of :math:`N` vertices which can
    belong to one of :math:`B` blocks, and the matrix :math:`e_{rs}` describes
    the number of edges from block :math:`r` to :math:`s` (or twice that number
    if :math:`r=s` and the graph is undirected).

    For an undirected graph, the number of distinct :math:`e_{rs}` matrices is
    given by,

    .. math::

       \Omega_m = \left(\!\!{B(B+1)/2 \choose E}\!\!\right)

    and for a directed graph,

    .. math::
       \Omega_m = \left(\!\!{B^2 \choose E}\!\!\right)


    where :math:`\left(\!{n \choose k}\!\right) = {n+k-1\choose k}` is the
    number of :math:`k` combinations with repetitions from a set of size
    :math:`n`. Hence, we have the description length of the edge counts

    .. math::

        -\ln P(\boldsymbol{e}) = \ln \Omega_m.

    For the node partition :math:`\boldsymbol{b}` we assume a two-level Bayesian
    hierarchy, where first the group size histogram is generated, and
    conditioned on it the partition, which leads to a description length:

    .. math::

       -\ln P(\boldsymbol{b}) = \ln\left(\!\!{B \choose N}\!\!\right) + \ln N! - \sum_r \ln n_r!,

    where :math:`n_r` is the number of nodes in block :math:`r`. If we forbid
    empty groups, with ``allow_empty == False``, this changes slightly to

    .. math::

       -\ln P(\boldsymbol{b}) = \ln {N - 1 \choose B - 1} + \ln N! - \sum_r \ln n_r!.

    The total information necessary to describe the model is then,

    .. math::

       -\ln P(\boldsymbol{e}, \boldsymbol{b}) = -\ln P(\boldsymbol{e}) - \ln P(\boldsymbol{b}).

    If ``nr`` is ``None``, it is assumed :math:`n_r=N/B`. If ``nr`` is
    ``False``, the partition term :math:`-\ln P(\boldsymbol{b})` is omitted
    entirely.

    References
    ----------

    .. [peixoto-parsimonious-2013] Tiago P. Peixoto, "Parsimonious module
       inference in large networks", Phys. Rev. Lett. 110, 148701 (2013),
       :doi:`10.1103/PhysRevLett.110.148701`, :arxiv:`1212.4794`.
    .. [peixoto-hierarchical-2014] Tiago P. Peixoto, "Hierarchical block
       structures and high-resolution model selection in large networks ",
       Phys. Rev. X 4, 011047 (2014), :doi:`10.1103/PhysRevX.4.011047`,
       :arxiv:`1310.4377`.
    .. [peixoto-model-2016] Tiago P. Peixoto, "Model selection and hypothesis
       testing for large-scale network models with overlapping groups",
       Phys. Rev. X 5, 011033 (2016), :doi:`10.1103/PhysRevX.5.011033`,
       :arxiv:`1409.3059`.

    """

    if directed:
        x = (B * B);
    else:
        x = (B * (B + 1)) / 2;
    if nr is False:
        L = lbinom(x + E - 1, E)
    else:
        L = lbinom(x + E - 1, E) + partition_entropy(B, N, nr, allow_empty)
    return L

def bethe_entropy(g, p):
    r"""Compute the Bethe entropy given the edge block membership marginals.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        The graph.
    p : :class:`~graph_tool.PropertyMap`
       Edge property map with edge marginals.

    Returns
    -------
    H : ``float``
        The Bethe entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_)
    Hmf : ``float``
        The "mean field" entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_),
        as would be returned by the :func:`mf_entropy` function.
    pv : :class:`~graph_tool.PropertyMap` (optional, default: ``None``)
        Vertex property map with vector-type values, storing the accumulated
        block membership counts. These are the node marginals, as would be
        returned by the :func:`collect_vertex_marginals` function.

    Notes
    -----

    The Bethe entropy is defined as,

    .. math::

        H = -\sum_{ij}A_{ij}\sum_{rs}\pi_{ij}(r,s)\ln\pi_{ij}(r,s) - \sum_i(1-k_i)\sum_r\pi_i(r)\ln\pi_i(r),

    where :math:`\pi_{ij}(r,s)` is the marginal probability that vertices
    :math:`i` and :math:`j` belong to blocks :math:`r` and :math:`s`,
    respectively, and :math:`\pi_i(r)` is the marginal probability that vertex
    :math:`i` belongs to block :math:`r`, and :math:`k_i` is the degree of
    vertex :math:`i` (or total degree for directed graphs).

    References
    ----------
    .. [mezard-information-2009] Marc Mézard, Andrea Montanari, "Information,
       Physics, and Computation", Oxford Univ Press, 2009.
       :DOI:`10.1093/acprof:oso/9780198570837.001.0001`
    """
    H = 0
    pv =  g.new_vertex_property("vector<double>")

    H, Hmf  = libinference.bethe_entropy(g._Graph__graph,
                                         _prop("e", g, p),
                                         _prop("v", g, pv))
    return H, Hmf, pv


def mf_entropy(g, p):
    r"""Compute the "mean field" entropy given the vertex block membership marginals.

    Parameters
    ----------
    g : :class:`~graph_tool.Graph`
        The graph.
    p : :class:`~graph_tool.PropertyMap`
        Vertex property map with vector-type values, storing the accumulated block
        membership counts.

    Returns
    -------
    Hmf : ``float``
        The "mean field" entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_).

    Notes
    -----

    The "mean field" entropy is defined as,

    .. math::

        H = - \sum_{i,r}\pi_i(r)ln\pi_i(r),

    where :math:`\pi_i(r)` is the marginal probability that vertex :math:`i`
    belongs to block :math:`r`.

    References
    ----------
    .. [mezard-information-2009] Marc Mézard, Andrea Montanari, "Information,
       Physics, and Computation", Oxford Univ Press, 2009.
       :DOI:`10.1093/acprof:oso/9780198570837.001.0001`
    """

    return libinference.mf_entropy(g._Graph__graph,
                                   _prop("v", g, p))
def microstate_entropy(h):
    r"""Compute microstate entropy given a histogram of partitions.

    Parameters
    ----------
    h : :class:`~graph_tool.inference.PartitionHist` (optional, default: ``None``)
        Partition histogram.

    Returns
    -------
    H : ``float``
        The microstate entropy value (in `nats <http://en.wikipedia.org/wiki/Nat_%28information%29>`_).

    Notes
    -----

    The microstate entropy is defined as,

    .. math::

        H = - \sum_{\boldsymbol b}p({\boldsymbol b})\ln p({\boldsymbol b}),

    where :math:`p({\boldsymbol b})` is observed frequency of partition
    :math:`{\boldsymbol b}`.

    References
    ----------
    .. [mezard-information-2009] Marc Mézard, Andrea Montanari, "Information,
       Physics, and Computation", Oxford Univ Press, 2009.
       :DOI:`10.1093/acprof:oso/9780198570837.001.0001`
    """

    return libinference.partitions_entropy(h)

from . overlap_blockmodel import *