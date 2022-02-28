import re

import numpy as np
import networkx as nx
from shapely.geometry.polygon import LinearRing

import drawSvg as draw

#some helper functions
def normal_length1(vector, side):
    """Normal vector which points usually outside of the molecule"""
    
    unit_vector = vector/np.linalg.norm(vector)
    unit_normal = np.array([-1*side*unit_vector[1], side*unit_vector[0]])
        
    return unit_normal
    
def calc_degree(vector1, vector2):
    """Formula for calculating cosA between to vectors"""
    
    cos_alpha = np.dot(vector1, vector2)/(np.linalg.norm(vector1)*np.linalg.norm(vector2))
    degree = np.arccos(cos_alpha)*180/np.pi
    return degree
    
def deg_to_rad(deg):
    return deg*np.pi/180

#all other main functions
def comprehend_string(string, orient):
    """
    
    string: kernel string
    orient: on default 1 which is clockwise orientation, but -1 is accepted as ccw

    Cuts up kernel string with a regex and saves the info into an array. Orientation is saved for each domain. By checking which bracket
    is paired with which, we save pairedness into a dictionary. A list of list is created which will hold the domain types and will be
    filled in later.
    """
    #separate domains/count
    pattern = re.compile('@?-?[0-9]* ?[a-z|A-Z|0-9|\*|\)]+[ \(]? ?\+?|\) ?\+?')
    dom_raw_list = pattern.findall(string)
    dom_count = len(dom_raw_list)
        
    #create empty data structures     
    struct_info_array = np.zeros((dom_count, 2), dtype = 'int') #3 cols rn, sb not needed here
    #empty angles indicated by -999
    for row in range(dom_count):
        struct_info_array[row, 0] = -999
        struct_info_array[row, 1] = orient
    name_list1 = []
    
    #loop through the raw domains, extract  
    for dom_i in range(len(dom_raw_list)):
        dom = dom_raw_list[dom_i]
        dom = dom.strip()
        
        #fill table
        if re.search("\s", dom):
            dom_parts = re.split("\s", dom)
            
            for segm in dom_parts:
                if segm[0] == "@": #col0 angle
                    segm = segm.lstrip("@")
                    struct_info_array[dom_i, 0] = int(segm) #col0 to write to
                    #struct_info_array[dom_i, 1] = np.sign(int(segm)) #col1 to write to
                    
                #elif segm is "+": #col2 strand break
                #    struct_info_array[dom_i, 2] = 1         #col2 to write to
                elif segm != "+": #else:
                    name_list1.append(segm) #only keep the name part further
        else:
            name_list1.append(dom)
        
    #pairedness dictionary
    paired_dict = {}
    
    for dom_i2 in range(len(dom_raw_list)):
        if re.search("\(", dom_raw_list[dom_i2]):
            bracket_count = 1
            for next_dom in range(dom_i2 + 1, len(dom_raw_list)):
                if re.search("\(", dom_raw_list[next_dom]):
                    bracket_count += 1
                elif re.search("\)", dom_raw_list[next_dom]):
                    bracket_count -= 1
                if bracket_count == 0:
                    paired_dict[dom_i2] = next_dom
                    break
    
    #name list of list / finding hairpin loops
    #paireds in dict, hairpin loop known, other are unpaired -> substructure from these?
    name_final_lol = []
    
    for dom_i3 in range(len(name_list1)):
        
        #paireds
        if re.search("\(", name_list1[dom_i3]):
            typ = "paired"
            name = name_list1[dom_i3].strip("(+ ")
        
        #correcting closing bracket name
        elif re.search("\)", name_list1[dom_i3]):
            typ = "paired"
            index_of_pair = list(paired_dict.keys())[list(paired_dict.values()).index(dom_i3)]
            if re.search("\*", name_list1[index_of_pair]):
                name = name_list1[index_of_pair].strip("\(").strip("\*")
            else:
                name = name_list1[index_of_pair].strip("\(") + "*"
        
        #all other are unpaireds
        else:
            typ = "unpaired"
            name = name_list1[dom_i3].strip("+ ")
        
        name_final_lol.append([name, typ])
        
    return dom_raw_list, struct_info_array, paired_dict, name_final_lol
	
	
def process_length_input(struct_info_array, length_input_dict):
    """
    
    struct_info_array: from comprehend_string
    length_input_dict: a dictionary, keys are domain indices, values are lengths

    Create a new column for struct_array, where lengths are going to be stored.
    """

    zeros_list = [0 for i in range(len(struct_info_array[:,1]))]
    
    #put into list format and append list to struct_info_array as new col
    for length in length_input_dict:
        zeros_list[length] = length_input_dict[length]
        
    zeros_list = np.array(zeros_list).reshape((len(zeros_list), 1))
    struct_info_array = np.append(struct_info_array, zeros_list, axis=1)
        
    return struct_info_array
	

def create_skeleton(dom_raw_list, paired_dict, name_final_lol):
    #append edges to graph
    G = nx.Graph()
    
    for dom_to_dom in range(len(name_final_lol)-1):
        if re.search("\+", dom_raw_list[dom_to_dom]) == None:
            G.add_edge(dom_to_dom, dom_to_dom+1, color='r') #domains that have a connection point are conn. with red edge
    for pair in paired_dict:
        G.add_edge(pair, paired_dict[pair], color='g') #paired domains are connected with green edge
    
    return G
	
	
def find_type_of_unpaired(G, paired_dict, name_final_lol):
    #index of ununpaireds - these don't need to be checked
    paired_node_list = []
    for node in paired_dict:
        paired_node_list.append(node)
        paired_node_list.append(paired_dict[node])
    
    #get unpaired too
    unpaired_node_list = []
    for node_index in range(len(name_final_lol)):
        if node_index not in paired_node_list:
            unpaired_node_list.append(node_index)
    
    #cycles of the skeleton graph
    cycle_list_G = nx.cycle_basis(G) #neat!
    multiloop_list = []
    
    #look through cycles, where 3cycles are hairpin, 5c. are bulge and longer ones nest multiloops
    for cycle in cycle_list_G:
        if len(cycle) == 3: #always hairpin
            for item in cycle:
                if item in unpaired_node_list:
                    name_final_lol[item][1] = "hairpin loop" #name altering
                    unpaired_node_list.remove(item) #trim unpaired list from known items
        elif len(cycle) > 3 & len(cycle) % 2 == 1: #current bulge def
            unpaired_here = []
            for item in cycle:
                if item in unpaired_node_list:
                    unpaired_here.append(item)
            if len(unpaired_here) == 1:
                name_final_lol[unpaired_here[0]][1] = "bulgeloop" #name altering
                unpaired_node_list.remove(unpaired_here[0]) #trim unpaired list from known items
         
        #find multiloops, the corresponding domains are their own list item in multiloop_list
        if len(cycle) >= 6:
            curr_multiloop = []
            for index in range(len(cycle)):
                if cycle[index] in unpaired_node_list:
                    
                    if cycle[index-3] in unpaired_node_list :
                        connections = [n for n in G.edges.data(nbunch=cycle[index-1])]
                        for edge in connections:
                            #only need to check paired edge between these2
                            if edge[1] == cycle[index-2] and edge[2]['color'] == 'g': 
                                name_final_lol[cycle[index-3]][1] = "multiloop"
                                name_final_lol[cycle[index]][1] = "multiloop" #altering names
                                if cycle[index] not in curr_multiloop:
                                    curr_multiloop.append(cycle[index])
                                if cycle[index-3] not in curr_multiloop:
                                    curr_multiloop.append(cycle[index-3])
            if len(curr_multiloop) > 0:
                multiloop_list.append(curr_multiloop)
                for dom in curr_multiloop:
                    unpaired_node_list.remove(dom) #trim unpaired list from known items
    
    #unhinged sequence
    for node_index in unpaired_node_list:
        #first neighbors
        direct_neigh = [n for n in G.neighbors(node_index)] #neighbors of current node
        
        #start from unhinged and go until all neighbors are unpaired
        if len(direct_neigh) == 1:
            name_final_lol[node_index][1] = "unhinged"
            
            next_neigh = direct_neigh
            before = [node_index]
            while len(next_neigh) == 1 and next_neigh[0] in unpaired_node_list:
                name_final_lol[next_neigh[0]][1] = "unhinged sequence"
                before.append(next_neigh[0])
                next_neigh = [n for n in G.neighbors(next_neigh[0])]
                next_neigh.remove(before[-2])
            for step in before:
                unpaired_node_list.remove(step) #bug with this line, cuts out 20 below, could implement it a bit diff
    
    #remaining is unknown but probably linear
    for remain in unpaired_node_list:
        name_final_lol[remain][1] = "ambiguous"
    
    return name_final_lol, multiloop_list, cycle_list_G

	
def resolve_domains_in_cycles(struct_info_array, paired_dict, skeleton_graph):
    """
    
    This function handles a lot of stuff, separated in inner functions. We take cycles and assume their geometry 
    with trigonometric rules and fill in the information from the string.
    """
    
    cycle_list_G = nx.cycle_basis(skeleton_graph)
    
    paired_node_list = []
    for node in paired_dict:
        paired_node_list.append(node)
        paired_node_list.append(paired_dict[node])
        
    def give_pair_of_domain(dom_index, paired_dict):
        if dom_index in paired_dict:
            return paired_dict[dom_index]
        elif dom_index in list(paired_dict.values()):
            index_of_pair = list(paired_dict.keys())[list(paired_dict.values()).index(dom_index)]
            return index_of_pair
        else:
            raise ValueError('not actually paired')
    
    def initiate_polygon_table(cycle, struct_info_array, paired_dict, paired_node_list):
        """
        
        Create table for said polygon, the angles for the polygon depend only on relative angles from the string
        """
        
        cycle.sort()
        polygon_array = np.zeros((int(polygon_sides), 3), dtype='float')
        paired_dist = 10 #could be changed
        
        curr_row = 0
        for cyc_index in range(len(cycle)):
            input_ang = struct_info_array[cycle[cyc_index], 0]  #rel. angle
            side = struct_info_array[cycle[cyc_index], 1]
            
            if cycle[cyc_index] not in paired_node_list: #at the start of an unpaired there is always a PA
                polygon_array[curr_row, 0] = cycle[cyc_index]   #remember dom name to write back to
                polygon_array[curr_row, 2] = struct_info_array[cycle[cyc_index], 2]  #dom length
                    
                curr_row += 1
                
                if cycle[cyc_index-1] not in paired_node_list: #unpaired-unpaired
                    if input_ang == -999:
                        poly_ang = 0
                    else:
                        poly_ang = 180 + input_ang*side
                    polygon_array[curr_row, 1] = poly_ang
                
                else:                                          #unpaired - paired(prev)
                    if input_ang == -999:
                        poly_ang = 0
                    else:
                        poly_ang = 90 + input_ang*side
                    polygon_array[curr_row, 1] = poly_ang
            
            else:
                pair_of_curr = give_pair_of_domain(cycle[cyc_index], paired_dict)
                
                if cycle[cyc_index - 1] not in paired_node_list: #paired-unpaired
                    if input_ang == -999:
                        poly_ang = 0
                    else:
                        poly_ang = 90 + input_ang*side
                    
                    polygon_array[curr_row, 0] = cycle[cyc_index]   #remember dom name to write back to
                    polygon_array[curr_row, 1] = poly_ang
                    polygon_array[curr_row, 2] = paired_dist  #dom length
                
                    curr_row += 1 #update row in polygon_array
                    
                elif cycle[cyc_index - 1] in paired_node_list and cycle[cyc_index - 1] != pair_of_curr: #paired-paired
                    if input_ang == -999:
                        poly_ang = 0
                    else:
                        poly_ang = np.abs(input_ang)
                    
                    polygon_array[curr_row, 0] = cycle[cyc_index]   #remember dom name to write back to
                    polygon_array[curr_row, 1] = poly_ang
                    polygon_array[curr_row, 2] = paired_dist  #dom length
                
                    curr_row += 1 #update row in polygon_array
        
        return polygon_array
    
    def indicate_polygon_case(polygon_array):
        """
        
        Based on arbitrary rules involving properties already in the polygon table, this function determines, 
        in which class the current loop belongs, hence with which algorithm we want to solve it.
        """
                
        angles_not_given = 0
        length_not_given = 0
        case1 = False
        for row in range(len(polygon_array[:,0])):
            if polygon_array[row-1, 2] == 0 and polygon_array[row-1, 1] == 0 and polygon_array[row, 1] == 0:
                case1 = True
            if polygon_array[row, 2] == 0:
                length_not_given += 1
            if polygon_array[row, 1] == 0:
                angles_not_given += 1

        #case1
        if case1 == True:
            return "case 1"
        else:
            if length_not_given >= 2:  #parallelness needs to be solved in case2calc function
                return "case 2"
            elif length_not_given == 1:
                if angles_not_given >= 2:
                    return "case 3"
                else:
                    return "case 1" #1l 0-1a
            else:  #all lengths known - like crossroads!
                return "case 4"
            
    def case1calc(polygon_array):
        """

        Calculates points of a vector polygon where at least 1 edge length is unknown, we start
        drawing from (0, 0) with the next domain after this choosen 'loose edge'. We assume 
        angles and length for every cell other than this loose end. Calculate the vertices 
        and calculate the last edge length and two last angles.
        """
        
        def resort_list_to_start_by_spec_value(dom_index, list_with_value):
            index = list_with_value.index(dom_index)
            sorted_list = []
            for i in range(len(list_with_value)):
                sorted_list.append(list_with_value[index-len(list_with_value)+i])
            return sorted_list
        
        n_poly = len(polygon_array[:,0])
        given_ang_list = []
        
        for row in range(n_poly):
            if polygon_array[row-1, 2] == 0 and polygon_array[row-1, 1] == 0 and polygon_array[row, 1] == 0:
                start = row
            if polygon_array[row, 1] != 0:
                given_ang_list.append(polygon_array[row, 1])
        
        #distribute remaining angles
        poly_angles_sum = (n_poly - 2) * 180
        
        remaining_average = (poly_angles_sum - sum(given_ang_list))/(n_poly - len(given_ang_list))
        
        for row in range(n_poly):
            if polygon_array[row, 1] == 0:
                polygon_array[row, 1] = remaining_average #set all unknown angles for this mean value
                
        #get length/angle mean ratio
        known_len_angle_sum = []
        for row in range(n_poly):
            if polygon_array[row-1, 2] != 0: #non-null length
                len_angle_sum_multip_here = (polygon_array[row-1, 1] + polygon_array[row, 1]) * polygon_array[row-1, 2]
                known_len_angle_sum.append(len_angle_sum_multip_here)
        len_angle_ratio = np.mean(known_len_angle_sum)
        
        #fill length col by using len_angle_ratio
        for row in range(n_poly):
            if polygon_array[row-1, 2] == 0: #null length! :)
                new_len = len_angle_ratio / (polygon_array[row-1, 1] + polygon_array[row, 1])
                polygon_array[row-1, 2] = new_len
        
        #get drawing order by resorting index list
        index_list = [n for n in range(n_poly)]
        resorted_list = resort_list_to_start_by_spec_value(start, index_list)
        resorted_list.pop() #remove last item which is the 'loose end' domain 
        
        #the table is filled and we can start calculating coords at domain 'start'
        coords = np.array([[0, 0]])
        
        for index_in_list in range(len(resorted_list)):
            ang_here = polygon_array[resorted_list[index_in_list], 1]
            len_here = polygon_array[resorted_list[index_in_list], 2]
            
            if index_in_list == 0:
                curr_coords = np.array([[len_here, 0]])
                abs_ang = 0
            else:
                abs_ang += 180-ang_here #abs angle is the outside angle of the polygon, adding them up 
                curr_coords = coords[-1] + np.array([[np.cos(deg_to_rad(abs_ang)) * len_here, 
                                                     np.sin(deg_to_rad(abs_ang)) * len_here]])
            
            coords = np.append(coords, curr_coords, axis=0)
        
        #the last prev_coords we get and [0, 0] adds the last side of the polygon
        len_last = np.linalg.norm(coords[-1])
        polygon_array[start-1, 2] = len_last
        
        #angle at start-1
        vector1 = coords[-2] - coords[-1]
        vector2 = coords[0] - coords[-1]
        angle = calc_degree(vector1, vector2)
        
        polygon_array[start-1, 1] = angle
        
        #angle at start
        vector1 = coords[-1] - coords[0]
        vector2 = coords[1] - coords[0]
        angle = calc_degree(vector1, vector2)
        
        polygon_array[start, 1] = angle
        
        return polygon_array
        
    def feed_polyarr_to_structarray(struct_info_array, polygon_array, cycle, paired_node_list):
        """Write back calculated relative angles and domain lengths(if appropriate) to SIA"""
        
        for poly_index in range(len(polygon_array[:,0])):
            curr_dom_index = int(polygon_array[poly_index, 0])
            poly_ang = polygon_array[poly_index, 1]
            poly_len = polygon_array[poly_index, 2]
            
            side = struct_info_array[curr_dom_index, 1]
            cyc_index = cycle.index(curr_dom_index)
            
            if curr_dom_index not in paired_node_list: #at the start of an unpaired there is always a polygon angle
                struct_info_array[curr_dom_index, 2] = poly_len #if unpaired, the length we write over
                
                if cycle[cyc_index-1] not in paired_node_list: #unpaired-unpaired
                    rel_angle = (poly_ang - 180) / side
                    struct_info_array[curr_dom_index, 0] = rel_angle
                
                else:                                          #unpaired - paired(prev)
                    rel_angle = (poly_ang - 90) / side
                    struct_info_array[curr_dom_index, 0] = rel_angle
            
            else:
                if cycle[cyc_index - 1] not in paired_node_list: #paired-unpaired
                    rel_angle = (poly_ang - 90) / side
                    struct_info_array[curr_dom_index, 0] = rel_angle
                    
                else: #it is always the right paired-paired, as we only put those in the table
                    rel_angle = poly_ang*side
                    struct_info_array[curr_dom_index, 0] = rel_angle
        
        return struct_info_array
    
    #MAIN LOOP
    for cycle in cycle_list_G:
        #treat hairpin and paired into paired separately ~ no polygon there
        polygon_sides = 0
        for item in cycle:
            if item in paired_node_list:
                polygon_sides += 0.5
            else:
                polygon_sides += 1
                
        if polygon_sides == 2: #special cases
            if len(cycle) == 3: #hairpin, nothing to do here
                pass
            elif len(cycle) == 4: #paired into paired, can fill in two 0 angles
                cycle.sort()      #the domains will be after this always the secomd and fourth in cycle
                struct_info_array[cycle[1], 0] = 0
                struct_info_array[cycle[3], 0] = 0

        elif polygon_sides >= 3: #has a polygon
            polygon_array = initiate_polygon_table(cycle, struct_info_array, paired_dict, paired_node_list)
            
            if indicate_polygon_case(polygon_array) == "case 1":
                polygon_array = case1calc(polygon_array)
            elif indicate_polygon_case(polygon_array) == "case 2":
                pass
            elif indicate_polygon_case(polygon_array) == "case 3":
                pass
            elif indicate_polygon_case(polygon_array) == "case 4":
                pass
            
            #add calculated info back to SIA
            struct_info_array = feed_polyarr_to_structarray(struct_info_array, polygon_array, cycle, paired_node_list)
            
    return struct_info_array
	
	
#determine order of traversal (order of domains to build up)
def determine_contstruct_order(skeleton_graph, paired_dict):
    cycle_list_G = nx.cycle_basis(skeleton_graph)
    
    paired_node_list = []
    for node in paired_dict:
        paired_node_list.append(node)
        paired_node_list.append(paired_dict[node])
    
    def find_cycle_where_index_belongs(dom_index, cycle_list_G): #now finds largest cycle
        found = []
        for cycle in cycle_list_G:
            if dom_index in cycle and len(cycle) > len(found):
                found = cycle
        return found
    
    def resort_list_to_start_by_spec_value(dom_index, list_with_value):
        index = list_with_value.index(dom_index)
        sorted_list = []
        for i in range(len(list_with_value)):
            sorted_list.append(list_with_value[index-len(list_with_value)+i])
        return sorted_list
    
    def give_pair_of_domain(dom_index, paired_dict):
        if dom_index in paired_dict:
            return paired_dict[dom_index]
        elif dom_index in list(paired_dict.values()):
            index_of_pair = list(paired_dict.keys())[list(paired_dict.values()).index(dom_index)]
            return index_of_pair
        else:
            raise ValueError('not actually paired')
    
    #
    traverse_order = [0]
    
    for node in traverse_order:
        current_in_travord = traverse_order.index(node)
        counter = 0 #helps with insert index
        
        #pair primary
        if node in paired_node_list:
            pair_first = give_pair_of_domain(node, paired_dict) #pair first!
            if pair_first not in traverse_order:
                traverse_order.insert(current_in_travord + 1, pair_first)
                counter += 1
        
        #cycle secondary
        its_cycle = find_cycle_where_index_belongs(node, cycle_list_G)
        if len(its_cycle) != 0: #if in cycle
            its_cycle.sort() #sort it first, bugfix
            #resort so current is at first place
            resorted_cycle = resort_list_to_start_by_spec_value(node, its_cycle)
            
            for item in resorted_cycle:
                if item not in traverse_order:
                    traverse_order.insert(current_in_travord + counter + 1, item)
                    counter += 1
                        
        #other neighbor tertiary
        neighbors = skeleton_graph.neighbors(node)
        for neighbor in neighbors:  #other neighbor third!
            if neighbor not in traverse_order:
                traverse_order.append(neighbor)

    return traverse_order
	

def stepwise_buildup(struct_info_array2, paired_dict, traverse_order, skeleton_graph):
    """
    
    Taking the more-or-less filled SIA, we go along the traversing/drawing order, assume lengths and angles for
    unpaired, out-of-loop domains and loop starters.
    """
    
    paired_dist = 10
    
    dom_count = len(traverse_order)
    coordinate_array = np.zeros((dom_count, 4))
    
    #add a new col for absolute angle values
    struct_info_array2 = np.append(struct_info_array2, [[0] for i in range(len(struct_info_array2[:,0]))], axis=1)
    
    #average over given lengths to get a default unpaired length
    n_given_len = 0
    sum_length = 0
    for length in struct_info_array2[:,2]:
        if length != 0:
            n_given_len += 1
            sum_length += length
    if n_given_len == 0:
        default_length = 30
    else:
        default_length = sum_length/n_given_len
    
    #paired list as before
    paired_node_list = []
    for node in paired_dict:
        paired_node_list.append(node)
        paired_node_list.append(paired_dict[node])
        
    def get_red_neigh(domain_index, skeleton_graph):
        connections = [n for n in skeleton_graph.edges.data(nbunch=domain_index)]
        red_neigh = []
        for edge in connections:
            if edge[2]['color'] is 'r':
                red_neigh.append(edge[1])
        return red_neigh
    
    def give_pair_of_domain(dom_index, paired_dict):
        if dom_index in paired_dict:
            return paired_dict[dom_index]
        elif dom_index in list(paired_dict.values()):
            index_of_pair = list(paired_dict.keys())[list(paired_dict.values()).index(dom_index)]
            return index_of_pair
        else:
            raise ValueError('not actually paired')
    
    def get_preferred_angle(domain_index, paired_node_list, skeleton_graph, traverse_order): #would be better with name_final
        neighbors = get_red_neigh(domain_index, skeleton_graph)
        if min(neighbors) == domain_index-1: #connected to prev dom
            if traverse_order.index(min(neighbors)) < traverse_order.index(domain_index):
                if min(neighbors) in paired_node_list and domain_index in paired_node_list:
                    pref_angle = 60

                elif min(neighbors) in paired_node_list and domain_index not in paired_node_list:
                    pref_angle = 90

                elif min(neighbors) not in paired_node_list and domain_index in paired_node_list:
                    pref_angle = 90

                elif min(neighbors) not in paired_node_list and domain_index not in paired_node_list:
                    pref_angle = 0

            elif traverse_order.index(min(neighbors)) > traverse_order.index(domain_index): #if backwards (SB case)
                if max(neighbors) in paired_node_list and domain_index in paired_node_list:
                    pref_angle = 60

                elif max(neighbors) in paired_node_list and domain_index not in paired_node_list:
                    pref_angle = 45

                elif max(neighbors) not in paired_node_list and domain_index in paired_node_list:
                    pref_angle = 45

                elif max(neighbors) not in paired_node_list and domain_index not in paired_node_list:
                    pref_angle = 180
            
        else:  #not connected to prev dom, so only to next dom
            if min(neighbors) in paired_node_list and domain_index in paired_node_list:
                pref_angle = 60

            elif min(neighbors) in paired_node_list and domain_index not in paired_node_list:
                pref_angle = 45

            elif min(neighbors) not in paired_node_list and domain_index in paired_node_list:
                pref_angle = 45

            elif min(neighbors) not in paired_node_list and domain_index not in paired_node_list:
                pref_angle = 180

        return pref_angle
    
    def calc_coords_of_paired2(domain_index, struct_info_array2, coordinate_array):
        side = struct_info_array2[domain_index, 1]
        x0p, y0p, x1p, y1p = coordinate_array[give_pair_of_domain(domain_index, paired_dict)]
        pair1vec = np.array([x1p - x0p, y1p - y0p])
        
        #direction of normal_vec has to be opposite of the default normal_vec, as that points outside and this inside
        x0, y0 = np.array([x1p, y1p]) + normal_length1(pair1vec, side) * paired_dist *-1  
        x1, y1 = np.array([x0, y0]) + pair1vec * -1
        
        return x0, y0, x1, y1
    
    def calculate_other_two_coords(x, y, abs_angle, length):
        abs_angle_rad = deg_to_rad(abs_angle)
        
        next_x = x + np.cos(abs_angle_rad) * length
        next_y = y + np.sin(abs_angle_rad) * length
        
        return next_x, next_y
    
    def take_opposite_angle(angle):
        if angle > 0:
            angle -= 180
        else:
            angle += 180
        return angle
    
    #MAIN LOOP   
    for domain_index in traverse_order:
        written_rows = traverse_order[0 : traverse_order.index(domain_index)] #rows before current in TO
        
        if domain_index == 0: #0th dom
            if struct_info_array2[0, 0] == -999:
                abs_angle = 0
            else:
                abs_angle = struct_info_array2[0, 0]
                struct_info_array2[0, 3] = abs_angle  #the absolute angle summation starts here
            
            if struct_info_array2[0, 2] == 0:
                length_here = default_length
            else:
                length_here = struct_info_array2[0, 2]
            
            x1,y1 = calculate_other_two_coords(0, 0, abs_angle, length_here)
            coordinate_array[0, 2] = x1
            coordinate_array[0, 3] = y1
            
        else: #all other doms
            #when drawing the second pair, just offset the first
            if domain_index in paired_node_list and give_pair_of_domain(domain_index, paired_dict) in written_rows:
                x0, y0, x1, y1 = calc_coords_of_paired2(domain_index, struct_info_array2, coordinate_array)
                coordinate_array[domain_index, 0] = x0
                coordinate_array[domain_index, 1] = y0
                coordinate_array[domain_index, 2] = x1
                coordinate_array[domain_index, 3] = y1
                
                abs_of_pair = struct_info_array2[give_pair_of_domain(domain_index, paired_dict), 3]
                abs_here = take_opposite_angle(abs_of_pair)
                struct_info_array2[domain_index, 3] = abs_here
            
            #others: unpaireds, first pairs
            else:
                side = struct_info_array2[domain_index, 1]
                neighbors = get_red_neigh(domain_index, skeleton_graph)
                
                #rel_angle cases
                if struct_info_array2[domain_index, 0] == -999:
                    prefer = get_preferred_angle(domain_index, paired_node_list, skeleton_graph, traverse_order)
                    rel_angle = side * prefer
                else:
                    rel_angle = struct_info_array2[domain_index, 0]
                    
                #dom length cases
                if struct_info_array2[domain_index, 2] == 0:
                    length_here = default_length
                else:
                    length_here = struct_info_array2[domain_index, 2]
                
                
                #for doms > 0 and previous dom didn't have strand break, x0, y0 is the same as x1, y1 of prev dom
                if len(neighbors) == 1:
                    if neighbors[0] < domain_index:
                        coordinate_array[domain_index, 0] = coordinate_array[neighbors[0], 2]
                        coordinate_array[domain_index, 1] = coordinate_array[neighbors[0], 3]
                        
                        abs_angle = struct_info_array2[neighbors[0], 3] + rel_angle
                        struct_info_array2[domain_index, 3] = abs_angle #fill out abs angle field

                        x0, y0 = coordinate_array[domain_index, 0], coordinate_array[domain_index, 1]
                        x1, y1 = calculate_other_two_coords(x0, y0, abs_angle, length_here)
                        coordinate_array[domain_index, 2] = x1
                        coordinate_array[domain_index, 3] = y1
                        
                        
                    elif neighbors[0] > domain_index:
                        coordinate_array[domain_index, 2] = coordinate_array[neighbors[0], 0]
                        coordinate_array[domain_index, 3] = coordinate_array[neighbors[0], 1]
                        
                        abs_angle = take_opposite_angle(struct_info_array2[neighbors[0], 3]) - rel_angle
                        struct_info_array2[domain_index, 3] = abs_angle #fill out abs angle field

                        x1, y1 = coordinate_array[domain_index, 2], coordinate_array[domain_index, 3] #other way around
                        x0, y0 = calculate_other_two_coords(x1, y1, abs_angle, length_here)
                        coordinate_array[domain_index, 0] = x0
                        coordinate_array[domain_index, 1] = y0
                        
                if len(neighbors) == 2:
                    for neigh in neighbors:
                        if neigh < domain_index:
                            if neigh in written_rows:
                                coordinate_array[domain_index, 0] = coordinate_array[domain_index - 1, 2]
                                coordinate_array[domain_index, 1] = coordinate_array[domain_index - 1, 3]

                        if neigh > domain_index:
                            if neigh in written_rows:
                                coordinate_array[domain_index, 2] = coordinate_array[domain_index + 1, 0]
                                coordinate_array[domain_index, 3] = coordinate_array[domain_index + 1, 1]

                    for neigh in neighbors:
                        if neigh < domain_index:
                            if neigh not in written_rows:
                                abs_angle = take_opposite_angle(struct_info_array2[domain_index + 1, 3]) - rel_angle #take the abs angle of the written neighbor
                                struct_info_array2[domain_index, 3] = abs_angle #fill out abs angle field

                                x1, y1 = coordinate_array[domain_index, 2], coordinate_array[domain_index, 3]
                                x0, y0 = calculate_other_two_coords(x1, y1, abs_angle, length_here)
                                coordinate_array[domain_index, 0] = x0
                                coordinate_array[domain_index, 1] = y0
                        if neigh > domain_index:
                            if neigh not in written_rows:
                                abs_angle = struct_info_array2[neighbors[0], 3] + rel_angle
                                struct_info_array2[domain_index, 3] = abs_angle #fill out abs angle field

                                x0, y0 = coordinate_array[domain_index, 0], coordinate_array[domain_index, 1]
                                x1, y1 = calculate_other_two_coords(x0, y0, abs_angle, length_here)
                                coordinate_array[domain_index, 2] = x1
                                coordinate_array[domain_index, 3] = y1
                                
    return coordinate_array
	

def create_color_list(name_final_lol, paired_dict, multiloop_list, palette):
    """
    
    color list for each domain where: pairs are same color, multiloops are same color, 
    neighbor domains never same color
    """
    if palette is "IBM":
        colors = ["#648FFF", "#785EF0", "#DC267F", "#FE6100", "#FFB000"]
        background_cols = ["#BACCFB", "#CAC1F3", "#DC9BBB", "#F1C2A6", "#F9E5B9"] #more 'faded' color for rectangle between pairs
    elif palette is "Wong":
        colors = ["#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
        background_cols = ["#BDBDBD", "#EAC267", "#A5CCE2", "#77A79A", "#E8E4A7", "#6595B1", "#D49E73", "#C7A8BA"]
    elif palette is "Magma":
        colors = ["#fcfdbf", "#fe9f6d", "#de4968", "#8c2981", "#3b0f70", "#000004"]
        background_cols = ["#FFFFE8", "#F7C9B1", "#DCA6B1", "#BBA0B8", "#AA93C5", "#A9A9A9"]
    elif palette is "Plasma":
        colors = ["#f0f921", "#fdb42f", "#ed7953", "#cc4778", "#9c179e", "#5c01a6", "#0d0887"]
        background_cols = ["#EDF194", "#F9DCA8", "#E8AE9B", "#CC91A6", "#AE7FAF", "#9C78B9", "#7E7CC1"]
    else:
        print("Palette not yet available")
        
    def give_pair_of_domain(dom_index, paired_dict):
        if dom_index in paired_dict:
            return paired_dict[dom_index]
        elif dom_index in list(paired_dict.values()):
            index_of_pair = list(paired_dict.keys())[list(paired_dict.values()).index(dom_index)]
            return index_of_pair
        else:
            raise ValueError('not actually paired')
    
    curr_color_i = 0
    for dom_index in range(len(name_final_lol)):
        color_i = curr_color_i % len(colors)
        
        if len(name_final_lol[dom_index]) == 2:
            name, domtype = name_final_lol[dom_index]

            if domtype is "paired":
                pair_index = give_pair_of_domain(dom_index, paired_dict)
                
                #index2 is line color
                name_final_lol[dom_index].append(colors[color_i])
                name_final_lol[pair_index].append(colors[color_i])
                
                #index3 is only for paired, the between pair background color
                name_final_lol[dom_index].append(background_cols[color_i])
                name_final_lol[pair_index].append(background_cols[color_i])
                
                curr_color_i += 1

            elif domtype is "multiloop":
                for loop in multiloop_list:
                    if dom_index in loop:
                        for loop_index in loop:
                            name_final_lol[loop_index].append(colors[color_i])
                            
            else:
                name_final_lol[dom_index].append(colors[color_i])
                
                curr_color_i += 1
    
    return name_final_lol
	
	
def lin_domain(x0, y0, x1, y1, side, name, color):
    
    #draw line from 5' to 3'
    p = draw.Line(x0, y0, x1, y1, stroke = color)
                  
    offset_len = 7
    text_size = 7
    
    #label has to be on its correct side
    vector = np.array([x1 - x0, y1 - y0])
    label_offset = normal_length1(vector, side) * offset_len
    x_text, y_text = np.array([x0,y0]) + vector/2 + label_offset
        
    t = draw.Text(name, text_size, x_text, y_text, fill = color)
    
    return p, t
	
	
def hairpin_loop(x0, y0, x1, y1, side, name, color):
    
    offset_len = 7
    text_size = 7
    
    #calculate Arc arguments
    r = np.sqrt((x1-x0)**2 + (y1-y0)**2) #rn radius is the distance between paired regions, will change
    
    #vector is pointing from point 0 -> point 1
    vector = np.array([x1 - x0, y1 - y0])
    vec_len = np.linalg.norm(vector)
    unit_normal = normal_length1(vector, side)
    
    #calculate center and label anchor with r and vector
    diag = np.sqrt(r**2 - (vec_len/2)**2)
    cx, cy = np.array([x0,y0])+ vector/2 + diag*unit_normal
    x_text, y_text = np.array([x0,y0])+ vector/2 + (diag+r+offset_len)*unit_normal
       
    vec_center_0 = np.array([x0-cx,y0-cy])
    vec_center_1 = np.array([x1-cx,y1-cy])
    
    #correction for 0 angle
    if np.sign(vec_center_0[1]) == 0:
        which_side0 = 1
    else:
        which_side0 = np.sign(vec_center_0[1])
    
    if np.sign(vec_center_1[1]) == 0:
        which_side1 = 1
    else:
        which_side1 = np.sign(vec_center_1[1])    
    
    #bit unclear about start/end degrees, but works
    if side == -1:
        startdeg = calc_degree(vec_center_0, np.array([1,0]))*which_side0
        enddeg = calc_degree(vec_center_1, np.array([1,0]))*which_side1
    elif side == 1:
        enddeg = calc_degree(vec_center_0, np.array([1,0]))*which_side0
        startdeg = calc_degree(vec_center_1, np.array([1,0]))*which_side1
        
    #draw
    p = draw.Arc(cx, cy, r, startdeg, enddeg,
        stroke=color, stroke_width=0.9, fill="none") #stroke_width will be also important
    
    #LABEL - if without name, no label is added
    t = draw.Text(name, text_size, x_text, y_text, fill=color) #label
        
    return p, t
	
	
def bulgeloop(x0, y0, x1, y1, side, name, color):
    
    offset_len = 7
    text_size = 7
    
    #calculate Arc arguments
    r = (np.sqrt((x1 - x0)**2 + (y1 - y0)**2))/2 #rn radius is the distance between paired regions, will change
    
    #vector is pointing from point 0 -> point 1
    vector = np.array([x1 - x0, y1 - y0])
    vec_len = np.linalg.norm(vector)
    unit_normal = normal_length1(vector, side)
    
    #calculate center and label anchor with r and vector
    cx, cy = np.array([x0, y0])+ vector/2
    x_text, y_text = np.array([x0, y0])+ vector/2 + (r + offset_len) * unit_normal
       
    vec_center_0 = np.array([x0 - cx, y0 - cy])
    vec_center_1 = np.array([x1 - cx, y1 - cy])
    
    #correction for 0 angle
    if np.sign(vec_center_0[1]) == 0:
        which_side0 = 1
    else:
        which_side0 = np.sign(vec_center_0[1])
    
    if np.sign(vec_center_1[1]) == 0:
        which_side1 = 1
    else:
        which_side1 = np.sign(vec_center_1[1])    
    
    #bit unclear about start/end degrees, but works
    if side == -1:
        startdeg = calc_degree(vec_center_0, np.array([1,0])) * which_side0
        enddeg = calc_degree(vec_center_1, np.array([1,0])) * which_side1
    elif side == 1:
        enddeg = calc_degree(vec_center_0, np.array([1,0])) * which_side0
        startdeg = calc_degree(vec_center_1, np.array([1,0])) * which_side1
        
    #draw
    p = draw.Arc(cx, cy, r, startdeg, enddeg,
        stroke=color, stroke_width=0.9, fill="none") #stroke_width will be also important
    
    #LABEL - if without name, no label is added
    t = draw.Text(name, text_size, x_text, y_text, fill=color) #label
        
    return p, t
	
	
def multiloop(d, ml, coordinate_array, struct_info_array2, name_final_lol):
    
    offset_len = 7
    text_size = 7
    
    side = struct_info_array2[ml[0], 1]
    colors = [name_final_lol[i][2] for i in ml]
    names = [name_final_lol[i][0] for i in ml]
    
    #taken from https://stackoverflow.com/a/47198877
    def find_center(p_ex1, p_ex2, centroid):
        x1, y1 = p_ex1
        x2, y2 = p_ex2
        x3, y3 = centroid
        dx, dy = x2 - x1, y2 - y1
        det = dx * dx + dy * dy
        a = (dy * (y3 - y1) + dx * (x3 - x1))/det
        
        return x1 + a * dx, y1 + a * dy
    
    #put points into shapely, get centroid
    coords = []
    
    for dom in ml:
        point0 = (coordinate_array[dom, 0], coordinate_array[dom, 1])
        point1 = (coordinate_array[dom, 2], coordinate_array[dom, 3])
        coords.append(point0)
        coords.append(point1)
        
    centroid = LinearRing(coords).centroid
    centroid = np.asarray(centroid)
    
    #get centerpoint for each arc
    for dom_i in range(len(ml)):
        dom = ml[dom_i]
        
        x0 = coordinate_array[dom, 0]
        y0 = coordinate_array[dom, 1]
        x1 = coordinate_array[dom, 2]
        y1 = coordinate_array[dom, 3]
        
        #vector is pointing from point 0 -> point 1
        vector = np.array([x1 - x0, y1 - y0])
        unit_normal = normal_length1(vector, side) #this should go to the other side as the hairpin center!
        
        p_ex1 = np.array([x0, y0]) + vector/2
        p_ex2 = np.array([x0, y0]) + vector/2 + unit_normal * -1
        
        cx, cy = find_center(p_ex1, p_ex2, centroid)
        radius = np.linalg.norm(np.array([x1 - cx, y1 - cy]))
        
        x_text, y_text = np.array([cx, cy]) + unit_normal * (radius + offset_len)
        
        vec_center_0 = np.array([x0 - cx, y0 - cy])
        vec_center_1 = np.array([x1 - cx, y1 - cy])

        #correction for 0 angle
        if np.sign(vec_center_0[1]) == 0:
            which_side0 = 1
        else:
            which_side0 = np.sign(vec_center_0[1])

        if np.sign(vec_center_1[1]) == 0:
            which_side1 = 1
        else:
            which_side1 = np.sign(vec_center_1[1])    

        #bit unclear about start/end degrees, but works
        if side == -1:
            startdeg = calc_degree(vec_center_0, np.array([1,0])) * which_side0
            enddeg = calc_degree(vec_center_1, np.array([1,0])) * which_side1
        elif side == 1:
            enddeg = calc_degree(vec_center_0, np.array([1,0])) * which_side0
            startdeg = calc_degree(vec_center_1, np.array([1,0])) * which_side1

        #append to draw object
        m = draw.Arc(cx, cy, radius, startdeg, enddeg,
            stroke=colors[dom_i], stroke_width=0.9, fill="none") #stroke_width will be also important

        #LABEL - if without name, no label is added
        t = draw.Text(names[dom_i], text_size, x_text, y_text, fill=colors[dom_i]) #label
        
        d.append(m)
        d.append(t)
    
    return d
	
	
def draw_image_from_coords(coordinate_array, struct_info_array2, name_final_lol, multiloop_list, paired_dict):
    
    def get_canvas_size(coordinate_array):
        x_coords = list(coordinate_array[:,0]) + list(coordinate_array[:,2])
        y_coords = list(coordinate_array[:,1]) + list(coordinate_array[:,3])
        
        x_len = (max(x_coords) - min(x_coords)) + 65 #canvas 14 incre. larger than max range
        y_len = (max(y_coords) - min(y_coords)) + 65
        
        origo = [min(x_coords) - 30, min(y_coords) - 30] #origin 7 increment lower - which is the text_offset
        canvas = [int(x_len), int(y_len)]
        
        return canvas, origo
    
    #draw process
    SCALE_METRIC = 3     #?
    CANVAS_SIZE, ORIGIN = get_canvas_size(coordinate_array)

    #image instancing
    d = draw.Drawing(CANVAS_SIZE[0], CANVAS_SIZE[1], origin=ORIGIN, displayInline=False)
    d.setPixelScale(SCALE_METRIC)
    
    for row in range(len(coordinate_array[:,0])):
        x0 = coordinate_array[row, 0]
        y0 = coordinate_array[row, 1]
        x1 = coordinate_array[row, 2]
        y1 = coordinate_array[row, 3]
        side = struct_info_array2[row, 1]
        name = name_final_lol[row][0]
        color = name_final_lol[row][2]
        
        #endings indicated
        if row == 0:
            dom1vec = np.array([x0 - x1, y0 - y1])
            vec5prime = dom1vec/np.linalg.norm(dom1vec) * 10
            t5 = draw.Text("5'", 6, vec5prime[0], vec5prime[1], fill='black')
            d.append(t5)
        elif row == len(coordinate_array[:,0])-1:
            dom_last_vec = np.array([x1 - x0, y1 - y0])
            vec3prime = np.array([x1,y1]) + dom_last_vec/np.linalg.norm(dom_last_vec) * 7
            t3 = draw.Text("3'", 6, vec3prime[0], vec3prime[1], fill='black')
            d.append(t3)
        
        #some dom types have special functions
        if name_final_lol[row][1] == "hairpin loop":
            p, t = hairpin_loop(x0, y0, x1, y1, side, name, color)
        elif name_final_lol[row][1] == "bulgeloop":
            p, t = bulgeloop(x0, y0, x1, y1, side, name, color)
        elif name_final_lol[row][1] != "multiloop":
            p, t = lin_domain(x0, y0, x1, y1, side, name, color)
        d.append(p)
        d.append(t)
    
    #multiloops done together
    for ml in multiloop_list:
        d = multiloop(d, ml, coordinate_array, struct_info_array2, name_final_lol)

    #rectangle between paired doms
    for pair1 in paired_dict:
        pair2 = paired_dict[pair1]
        color = name_final_lol[pair1][3]
        x0, y0, x1, y1 = coordinate_array[pair1]
        x0p, y0p, x1p, y1p = coordinate_array[pair2]
        
        r = draw.Lines(x0, y0, x1, y1,
                   x0p, y0p, x1p, y1p,
                   close = False, fill = color)
        d.append(r)

    return d
	
	
def domain_visualization(string, length_input_dict={}, palette="IBM", filename=None, orient=1):
    """
    
    Includes all previously defined functions, save intermediate data structures inside as local variables and 
    renders the drawing.
    """
    
    dom_raw_list, struct_info_array, paired_dict, name_final_lol = comprehend_string(string, orient)
    struct_info_array = process_length_input(struct_info_array, length_input_dict)
    skeleton_graph = create_skeleton(dom_raw_list, paired_dict, name_final_lol)
    name_final_lol, multiloop_list, cycle_list_G = find_type_of_unpaired(skeleton_graph, paired_dict, name_final_lol)
    
    #if no cycles skip cycle resolve:
    if len(cycle_list_G) == 0:
        struct_info_array2 = struct_info_array
    else:
        struct_info_array2 = resolve_domains_in_cycles(struct_info_array, paired_dict, skeleton_graph)
    
    traverse_order = determine_contstruct_order(skeleton_graph, paired_dict)
    coordinate_array = stepwise_buildup(struct_info_array2, paired_dict, traverse_order, skeleton_graph)    
    
    name_final_lol = create_color_list(name_final_lol, paired_dict, multiloop_list, palette)
    d = draw_image_from_coords(coordinate_array, struct_info_array2, name_final_lol, multiloop_list, paired_dict)
    
    if filename != None:
        d.saveSvg(filename)
	
    return d