module elements_nodes_neighbours

  use data_structure

  type (type_node_list)        :: node_list
  type (type_node_list),pointer:: aux_node_list
  type (type_element_list)     :: element_list
  integer,allocatable          :: element_neighbours(:,:)

end module
