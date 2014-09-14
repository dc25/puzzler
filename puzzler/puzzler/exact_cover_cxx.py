#!/usr/bin/env python
# $Id$

# Author: David Goodger <goodger@python.org>
# Copyright: (C) 1998-2012 by David J. Goodger
# License: GPL 2 (see __init__.py)

"""
A wrapper around the 'exactcover_cxx' C++ extension by David Compton

The 'exactcover_cxx' C++ extension implements Donald E. Knuth's 
'Algorithm X' [1]_ for the generalized exact cover problem [2]_ 
using the 'Dancing Links' technique [3]_ ('DLX').

.. [1] http://en.wikipedia.org/wiki/Knuth%27s_Algorithm_X
.. [2] http://en.wikipedia.org/wiki/Exact_cover
.. [3] http://en.wikipedia.org/wiki/Dancing_Links
"""

import exactcover_cxx
from puzzler.utils import thousands


class ExactCover(object):

    """
    Given a sparse matrix of 0s and 1s, find every set of rows containing
    exactly one 1 in each primary column (and at most one 1 in each secondary
    column).  See `load_matrix` for a description of the data structure.
    Uses the 'exactcover_cxx' C++ extension implementation of the Dancing Links
    approach to Knuth's Algorithm X.
    """

    def __init__(self, matrix=None, secondary=0, state=None):
        """
        Parameters:

        * `matrix` & `secondary`: see `self.load_matrix`.

        * `state`: a `puzzler.SessionState` object which stores the runtime
          state of this puzzle (we're resuming a previously interrupted
          puzzle), or None (no state, we're starting from the beginning).
        """
        self.solver = None
        """An `exactcover_cxx.Coverings` object, set in
        `self.load_matrix()`."""

        self.solution = []
        self.num_solutions = 0
        self._num_previous_searches = 0

        """Keeps track of the total number of searches in previous
        sub-puzzles, since the cxx extension doesn't."""

        if state:
            self.solution = state.solution
            self.num_solutions = state.num_solutions 
            self._num_previous_searches = state.num_searches

        if matrix:
            self.load_matrix(matrix, secondary)

    def load_matrix(self, matrix, secondary=0):
        """
        Convert the input `matrix` into a form compatible with the
        `exactcover_cxx` C++ extension and load it into an 
        `exactcover_cxx.Coverings` object.

        The input `matrix` is a two-dimensional list of lists:

        * Each row is a list of equal length.

        * The first row contains the column names: first the puzzle piece
          names, then the solution space coordinates.  For example::

              ['A', 'B', 'C', '0,0', '1,0', '0,1', '1,1']

        * The subsequent rows consist of 1 & 0 (True & False) values.  Each
          row contains a 1/True value in the column identifying the piece, and
          1/True values in each column identifying the position.  There must
          be one row for each possible position of each puzzle piece.

        The `secondary` parameter is the number of secondary (rightmost)
        columns: columns which may, but need not, participate in the solution.

        The converted data structure consists of a list of lists of column
        names.
        """
        if self.solver:
            self._num_previous_searches += self.solver.num_searches
        rows = [ [ item for item in row if item ] for row in matrix[1:]]

        self.solver = exactcover_cxx.Coverings( rows, matrix[0], self.solution, secondary)

    def getState(self):
        """
        Query the solver to get the current state.
        """
        state = self.solver.getState()
        return [state.getRow(r) for r in range(state.size())]

    def full_solution(self):
        """
        Required by report_search_state() in __init__.py
        """
        return self.solution

    def solve(self, level=0):
        """
        A generator that produces all solutions via the 
        `exactcover_cxx.Coverings` solver, using Algorithm X.
        """
        while( True ):
            sol = self.solver.getSolution()
            if sol is None:
                # This state gets passed to solver if solver is restarted.
                self.solution = self.getState()
                return;
            self.solution = [sol.getRow(r) for r in range(sol.size())]
            yield self.solution

    def format_solution(self):
        """Return a simple formatted string representation of the solution."""
        self.num_solutions += 1
        parts = [
            'solution %i (%s searches):'
            % (self.num_solutions, thousands(self.num_searches))]
        for row in self.solution:
            parts.append(
                ' '.join(cell for cell in row
                         # omit secondary columns (intersections):
                         if not ((',' in cell) and (cell.endswith('i')))))
        return '\n'.join(parts)

    @property
    def num_searches(self):
        """The number of search operations tried so far."""
        if self.solver:
            return self._num_previous_searches + self.solver.num_searches
        else:
            return self._num_previous_searches
