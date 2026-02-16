module c_interface
  !! Fortran to C interface for building dynamically sized edge connectivity arrays.

  interface
     subroutine create_adj_list(np, nedges, nedno1, nedno2) bind(C, name='create_adj_list')
  	   !! Create an adjacency list for each node, i.e. a list of its
  	   !! neighbors.
  	   !!
  	   !! The structure is a list of lists of pairs. Each node gets
  	   !! a list of neighbors. Each neighbor is represented with a pair: the
  	   !! adjacent node index and the edge index. Note that these are converted
  	   !! to C's 0-first convention under-the-hood.
       use, intrinsic :: iso_c_binding, only : c_int

       implicit none

       integer(c_int), value :: np
	   !! Number of nodes in the mesh.
	   integer(c_int), value :: nedges
	   !! Number of edges in the mesh.
       integer(c_int) :: nedno1(nedges)
	   !! List of starting node indexes for each edge.
	   integer(c_int) :: nedno2(nedges)
	   !! List of terminal node indexes for each edge.

     end subroutine create_adj_list

     integer function get_edge_no(n1, n2) bind(C, name='get_edge_no')
  	   !! Look up the edge index between a node and its neighbor.
  	   !!
  	   !! The Fortran 1-first convention should be used when providing
  	   !! arguments. Conversion to C's convention is done under-the-hood.
  	   !!
  	   !! @warning
  	   !! If nothing is found, 0 is returned.
  	   !! @endwarning
       use, intrinsic :: iso_c_binding, only : c_int
       implicit none

       integer(c_int), value :: n1
	   !! The first node index.
	   integer(c_int), value :: n2
	   !! The second node index.

     end function get_edge_no

  end interface
end module c_interface
