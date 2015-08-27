
#=======================================================================================================================
# Coarse detection of the worm posture.
#=======================================================================================================================

#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: initializedcheck=False


cimport cython
from cython cimport view
from cython.operator cimport dereference as deref, preincrement as inc
from libcpp cimport bool
from libc.math cimport log
cdef extern from "math.h":
    double INFINITY
from libcpp.pair cimport pair
from libcpp.vector cimport vector
from libcpp.map cimport map
from libcpp.set cimport set as cset

import numpy as np
cimport numpy as cnp
import scipy.ndimage

##import plot



cdef pair[int,int] *DIRECTIONS = [(0,1), (1,1), (1,0), (1,-1), (0,-1), (-1,-1), (-1,0), (-1,1)]


pyORIENTATIONS = {(0,1): 0, (1,1): 1, (1,0): 2, (1,-1): 3, (0,-1): 4, (-1,-1): 5, (-1,0): 6, (-1,1): 7}
cdef map[pair[int,int], int] ORIENTATIONS = pyORIENTATIONS


cdef struct detection_parameters:
    cnp.uint8_t[:,:,::1] features
    double log_pobj, log_1mpobj, log_pneighbor, log_1mpneighbor, log_pbg, log_1mpbg, lam, A, tau
    int search_depth, search_breadth
    bool remote

cdef detection_parameters DETECTION_PARAMETERS


 
cdef int _intersection(vector[pair[int,int]] & v1, vector[pair[int,int]] & v2):
    """
    Computes the number of elements in the intersection of two vectors (of integer pairs).
    
    Parameters
    ----------
    v1, v2 : vectors of integer pairs
        Input vectors, supposed to have no duplicates.
        
    Returns
    -------
    out : int
        Number of elements in the intersection. 
    """

    # create set from first vector
    cdef cset[pair[int,int]] s1
    cdef vector[pair[int,int]].iterator it1 = v1.begin()
    while it1 != v1.end():
        s1.insert(deref(it1))
        inc(it1)
    
    # loop over elements in second vector and count how often the element is in the set 
    cdef int res = 0
    cdef vector[pair[int,int]].iterator it2 = v1.begin()
    while it2 != v1.end():
        res += s1.count(deref(it2))
        inc(it2)
    
    return res


cdef _set_detection_parameters(cnp.ndarray[dtype=cnp.uint8_t,ndim=3,mode='c'] features, \
                               double lam, double A, \
                               double pobj, double pneighbor, double pbg, \
                               int search_depth, int search_breadth, \
                               bool remote, double tau=-1):
    """
    Stores the detection parameters for global access.
    
    Parameters
    ----------
    features : array
        Detected worm features.
    lam : double
        Expected length of posture.
    A : double
        Penalty for deviation from expected length.
    pobj, pneighbor, pbg : double
        Probability of feature on object/neighbor/background.
    search_depth : int
        Nr of search steps w/o improvement before stop.
    search_breadth : int
        Size of candidate list before remote search is turned off.
    remote : bool
        Whether an extended search should be used for remote candidates.
    tau : bool, optional
        Threshold for log-likelihood ratio before candidates are no longer followed. If -1, 2 times
        the maximal contrast (log-likelihood for optimal step minus log-likelihood for worst step) is used.
    """
    
    # create one-point border around features to avoid getting out of bounds
    DETECTION_PARAMETERS.features = np.zeros((features.shape[0]+2, features.shape[1]+2, 4), dtype=np.uint8)
    cdef int r, c, e 
    for r in range(features.shape[0]):
        for c in range(features.shape[1]):
            for e in range(4):
                DETECTION_PARAMETERS.features[r+1, c+1, e] = features[r,c,e]
    
    DETECTION_PARAMETERS.log_pobj = log(pobj)
    DETECTION_PARAMETERS.log_1mpobj = log(1 - pobj)
    DETECTION_PARAMETERS.log_pneighbor = log(pneighbor)
    DETECTION_PARAMETERS.log_1mpneighbor = log(1 - pneighbor)
    DETECTION_PARAMETERS.log_pbg = log(pbg)
    DETECTION_PARAMETERS.log_1mpbg = log(1 - pbg)
    DETECTION_PARAMETERS.lam = lam
    DETECTION_PARAMETERS.A = A
    
    if tau == -1:
        tau = (log(pobj)-log(pbg))+2*(log(pneighbor)-log(pbg))-(log(1-pobj)-log(1-pbg))-2*(log(1-pneighbor)-log(1-pbg))
        tau = 2 * tau
    
    DETECTION_PARAMETERS.tau = tau
    DETECTION_PARAMETERS.search_depth = search_depth
    DETECTION_PARAMETERS.search_breadth = search_breadth
    DETECTION_PARAMETERS.remote = remote


cdef double _log_likelihood_change(pair[int,int] last_point, int orientation, cnp.uint8_t[:,::1] blocked, \
                                   cnp.uint8_t[:,::1] other_blocked=None):
        """
        Computes the change in log-likelihood when the specified move is performed.
        
        Parameters
        ----------
        last_point : pair[int,int]
            Starting point of the move.
        orientation : int
            Orientation of the move.
        blocked : array
            Points blocked by the own subinstantiation. Not supposed to be None.
        other_blocked : array, optional
            Points blocked by the subinstantiation in the other direction. Can be None.
            
        Returns
        -------
        loglik_change : double
            Change in log-likelihood.
        """
        
        cdef pair[int,int] direction = DIRECTIONS[orientation]
        cdef pair[int,int] new_point = pair[int,int](last_point.first+direction.first, \
                                                     last_point.second+direction.second)
        
        if blocked[new_point.first, new_point.second] \
        or (other_blocked is not None and other_blocked[new_point.first, new_point.second]):
            return -INFINITY
            '''if new_point.first == 0 or new_point.first == DETECTION_PARAMETERS.features.shape[0]-1 \
            or new_point.second == 0 or new_point.second == DETECTION_PARAMETERS.features.shape[1]-1:
                return -INFINITY
            return 0'''
    
        cdef int orth_orientation = (orientation + 2) % 4
        cdef pair[int,int] orth_direction = DIRECTIONS[orth_orientation]
    
        # likelihood-ratio at new point
        cdef double loglik_change
        if DETECTION_PARAMETERS.features[new_point.first, new_point.second, orth_orientation]:
            loglik_change = DETECTION_PARAMETERS.log_pobj - DETECTION_PARAMETERS.log_pbg
        else:
            loglik_change = DETECTION_PARAMETERS.log_1mpobj - DETECTION_PARAMETERS.log_1mpbg

        # likelihood-ratio at neighbors (blocked neighbors do not contribute)
        cdef pair[int,int] neighbor
        if orth_orientation % 2 == 0: # vertical or horizontal move
            neighbor = pair[int,int](new_point.first+orth_direction.first, new_point.second+orth_direction.second)
            if not blocked[neighbor.first, neighbor.second] \
            and (other_blocked is None or not other_blocked[neighbor.first, neighbor.second]):
                if DETECTION_PARAMETERS.features[neighbor.first, neighbor.second, orth_orientation]:
                    loglik_change += DETECTION_PARAMETERS.log_pneighbor - DETECTION_PARAMETERS.log_pbg
                else:
                    loglik_change += DETECTION_PARAMETERS.log_1mpneighbor - DETECTION_PARAMETERS.log_1mpbg
            neighbor = pair[int,int](new_point.first-orth_direction.first, new_point.second-orth_direction.second)
            if not blocked[neighbor.first, neighbor.second] \
            and (other_blocked is None or not other_blocked[neighbor.first, neighbor.second]):
                if DETECTION_PARAMETERS.features[neighbor.first, neighbor.second, orth_orientation]:
                    loglik_change += DETECTION_PARAMETERS.log_pneighbor - DETECTION_PARAMETERS.log_pbg
                else:
                    loglik_change += DETECTION_PARAMETERS.log_1mpneighbor - DETECTION_PARAMETERS.log_1mpbg
        else: # diagonal move
            neighbor = pair[int,int](last_point.first+direction.first, last_point.second)
            if not blocked[neighbor.first, neighbor.second] \
            and (other_blocked is None or not other_blocked[neighbor.first, neighbor.second]):
                if DETECTION_PARAMETERS.features[neighbor.first, neighbor.second, orth_orientation]:
                    loglik_change += DETECTION_PARAMETERS.log_pneighbor - DETECTION_PARAMETERS.log_pbg
                else:
                    loglik_change += DETECTION_PARAMETERS.log_1mpneighbor - DETECTION_PARAMETERS.log_1mpbg
            neighbor = pair[int,int](last_point.first, last_point.second+direction.second)
            if not blocked[neighbor.first, neighbor.second] \
            and (other_blocked is None or not other_blocked[neighbor.first, neighbor.second]):
                if DETECTION_PARAMETERS.features[neighbor.first, neighbor.second, orth_orientation]:
                    loglik_change += DETECTION_PARAMETERS.log_pneighbor - DETECTION_PARAMETERS.log_pbg
                else:
                    loglik_change += DETECTION_PARAMETERS.log_1mpneighbor - DETECTION_PARAMETERS.log_1mpbg
        
        return loglik_change


cdef class Instantiation:
    """
    Represents posture instantiations.
    """
    
    cdef:
        double log_likelihood
        double log_posterior
        vector[pair[int,int]] points
        cnp.uint8_t[:,::1] blocked
        int last_orientation
        int last_turn
    
    
    def __init__(self, Instantiation instantiation=None):
        """
        Initialize instantiation.
        
        Parameters
        ----------
        instantiation : Instantiation, optional
            Instantiation to copy.
        """
        
        if instantiation is not None:
            self.log_likelihood = instantiation.log_likelihood
            self.log_posterior = instantiation.log_posterior
            self.points = instantiation.points
            self.blocked = view.array(shape=(DETECTION_PARAMETERS.features.shape[0],\
                                             DETECTION_PARAMETERS.features.shape[1]), \
                                      itemsize=sizeof(cnp.uint8_t), format="c")
            self.blocked[:] = instantiation.blocked
            #self.blocked = np.copy(instantiation.blocked)
            self.last_orientation = instantiation.last_orientation
            self.last_turn = instantiation.last_turn
        else:
            self.log_likelihood = -INFINITY
            self.log_posterior = -INFINITY
            self.blocked = None
            self.last_orientation = -1


    cdef Instantiation _make_all_admissible_moves1(self, candidates, Instantiation best_instantiation, \
                                                   cnp.uint8_t[:,::1] other_blocked=None):
        """
        Tries all admissible moves (based on the last move and blocked points) and updates the candidate dictionary
        as well as the best instantiation.
        
        Parameters
        ----------
        candidates : dictionary
            Dictionary of candidates, to be updated.
        best_instantiation : Instantiation
            The currently best instantiation.
        other_blocked : array, optional
            Points blocked by the subinstantiation in the other direction.
            
        Returns
        -------
        best_instantiation : Instantiation
            The new best instantiation.
        """
        
        cdef int orientation, i
        cdef pair[int,int] direction, orth_direction
        cdef pair[int,int] last_point = self.points.back(), new_point
        cdef double loglik_change, new_loglik, new_logpost
    
        cdef cnp.uint8_t[:,::1] new_blocked
        cdef vector[pair[int,int]] new_points
        cdef Instantiation candidate
    
        # try all possible moves (left, straight, right)
        cdef int turn
        for turn in range(-1,1+1):
            
            # no consecutive turns in the same direction
            if turn == -1 and self.last_turn == -1: continue
            if turn == 1 and self.last_turn == 1: continue
    
            # new log-likelihood
            orientation = (self.last_orientation + turn + 8) % 8  # avoids negative modulus
            direction = DIRECTIONS[orientation]
            loglik_change = _log_likelihood_change(last_point, orientation, self.blocked, other_blocked)
            if loglik_change == -INFINITY:
                continue
            new_loglik = self.log_likelihood + loglik_change
            
            # store resulting instantiation as candidate
            # (if the new point was not already reached with a better likelihood)
            py_new_point = (last_point.first+direction.first, last_point.second+direction.second)
            if py_new_point not in candidates \
                or new_loglik > (<Instantiation>(candidates[py_new_point])).log_likelihood:
      
                new_point = pair[int,int](last_point.first+direction.first, last_point.second+direction.second)
      
                # update candidate
                candidate = Instantiation()
                candidate.log_likelihood = new_loglik
                candidate.points = vector[pair[int,int]](self.points)
                candidate.points.push_back(new_point)
                candidate.blocked = view.array(shape=(DETECTION_PARAMETERS.features.shape[0],\
                                                      DETECTION_PARAMETERS.features.shape[1]), \
                                               itemsize=sizeof(cnp.uint8_t), format="c")
                candidate.blocked[:] = self.blocked
                for i in range(2, 7):  # block all neighbors except straight, left and right turn
                    candidate.blocked[last_point.first + DIRECTIONS[(orientation+i)%8].first, \
                                      last_point.second + DIRECTIONS[(orientation+i)%8].second] = True
                candidate.blocked[new_point.first, new_point.second] = True
                candidate.last_orientation = orientation
                candidate.last_turn = turn
                
                # add candidate to list
                candidates[py_new_point] = candidate
                
                # update best instantiation
                if new_loglik > best_instantiation.log_likelihood:
                    best_instantiation = candidates[py_new_point]
                    
        return best_instantiation

    
    cdef update_posterior(self):
        self.log_posterior = self.log_likelihood - \
                             DETECTION_PARAMETERS.A * (self.points.size() - DETECTION_PARAMETERS.lam)**2
    

cdef _detect_posture1(Instantiation init_instantiation, bool doviz=False):
    """
    Runs a posture detection starting from the provided instantiation (using the current detection setting).
    
    Parameters
    ----------
    init_instantiation : Instantiation
        Instantiation to start from.
    doviz : bool, optional
        Whether a visualization should be created.
        
    Returns
    -------
    best_instantiation : Instantiation
        The detected posture.
    viz : array
        Visualizations after 15, 30, 45 iterations, and the final detection.
    """

    # initialize best instantiation
    cdef Instantiation best_instantiation = init_instantiation

    # initialize candidate list
    candidates = {init_instantiation.points.back(): best_instantiation}
    
    ### run detection (until no improvement happened for a while) ###
    
    cdef int length = init_instantiation.points.size()
    cdef int last_improvement = length
    cdef double previous_length_best_loglik = init_instantiation.log_likelihood
    cdef Instantiation candidate, current_length_best_instantiation
    if doviz:
        visited = set(candidates.keys())
        viz = np.zeros((4,DETECTION_PARAMETERS.features.shape[0],DETECTION_PARAMETERS.features.shape[1],3), \
               dtype=np.uint8)
    else:
        viz = None
    
    while True:
        
        # setup for next length
        length += 1
        new_candidates = {}
        current_length_best_instantiation = Instantiation()

        # loop through candidates
        for candidate in candidates.values():
            
            # if above the log-likelihood threshold...
            # ...or candidate is remote and above lower threshold (and candidate list not too long)
            if candidate.log_likelihood > previous_length_best_loglik - DETECTION_PARAMETERS.tau \
                or (DETECTION_PARAMETERS.remote \
                    and candidate.log_likelihood > previous_length_best_loglik - 2*DETECTION_PARAMETERS.tau \
                    and len(candidates) <= DETECTION_PARAMETERS.search_breadth) \
                    and length - _intersection(best_instantiation.points,candidate.points) \
                        >= DETECTION_PARAMETERS.lam//5:
                # make all moves and update the new best instantiation                
                current_length_best_instantiation = \
                    candidate._make_all_admissible_moves1(new_candidates, current_length_best_instantiation)
        
        # check if the best instantiation of this iteration has better posterior
        if current_length_best_instantiation.log_likelihood - \
        DETECTION_PARAMETERS.A * (length - DETECTION_PARAMETERS.lam)**2 > best_instantiation.log_posterior:
            best_instantiation = current_length_best_instantiation
            best_instantiation.log_posterior = current_length_best_instantiation.log_likelihood - \
                                               DETECTION_PARAMETERS.A * (length - DETECTION_PARAMETERS.lam)**2 
            last_improvement = length
        
        # store results from current iteration
        previous_length_best_loglik = current_length_best_instantiation.log_likelihood
        candidates = new_candidates
        if doviz:
            visited = visited.union(set(new_candidates.keys()))
            viz_image = -1
            if length == 15:
                viz_image = 0
            elif length == 30:
                viz_image = 1
            elif length == 45:
                viz_image = 2
            if viz_image != -1:
                rr, cc = zip(*visited)
                viz[viz_image, rr, cc] = (0,128,128)
                rr, cc = np.where(np.any(DETECTION_PARAMETERS.features,2))
                viz[viz_image, rr,cc,] = 128
                rr, cc = zip(*candidates.keys())
                viz[viz_image, rr, cc] = (128,0,128)
                rr, cc = zip(*best_instantiation.points)
                viz[viz_image, rr,cc,] = 255             
        
        # stop search if no improvement has happened for too long
        if length == last_improvement + DETECTION_PARAMETERS.search_depth:
            break
    
    if doviz:
        rr, cc = zip(*visited)
        viz[3, rr, cc] = (0,128,128)
        rr, cc = np.where(np.any(DETECTION_PARAMETERS.features,2))
        viz[3, rr,cc,] = 128
        rr, cc = zip(*best_instantiation.points)
        viz[3, rr,cc,] = 255
        
    return best_instantiation, viz


cdef _find_border_points(A):
    """
    Finds a point on the top, bottom, left and right border of the binary array.
    
    Parameters
    ----------
    A : array
        Binary input array.
        
    Returns
    -------
    starting_points : list of integer pairs
        The four border points; or None if A is zero.
    """
    
    rs, cs = np.nonzero(A)
    if rs.size == 0:
        return None
    
    # top
    r_top = np.min(rs)
    c_top = np.int(np.round(np.mean(cs[rs==r_top])))
    top = (r_top, c_top)
    
    # bottom
    r_bottom = np.max(rs)
    c_bottom = np.int(np.round(np.mean(cs[rs==r_bottom])))
    bottom = (r_bottom, c_bottom)
    
    # left
    c_left = np.min(cs)
    r_left = np.int(np.round(np.mean(rs[cs==c_left])))
    left = (r_left, c_left)
    
    # bottom
    c_right = np.max(cs)
    r_right = np.int(np.round(np.mean(rs[cs==c_right])))
    right = (r_right, c_right)
    
    return [top, bottom, left, right]


def detect_worm1(cnp.ndarray[cnp.uint8_t, ndim=3, mode='c'] head_features, double head_lam, double head_A, \
                 cnp.ndarray[cnp.uint8_t, ndim=3, mode='c'] worm_features, double worm_lam, double worm_A, \
                 double pobj = .7, double pneighbor = .5, double pbg = .3, \
                 int soft_min_size=10, int hard_min_size=5, double min_logpost=90, cnp.uint8_t verbose=False):
    """
    First detects the head by running posture detections from reasonable start points
    (provided that the connected components of the head features are valid).
    Then detects the worm body by starting separate searches from both ends of the detected head.
    
    Parameters
    ----------
    head_features, worm_features : array
        Detected head/worm features.
    head_lam, worm_lam : double
        Expected length of head/worm posture (in coarse blocks). 
    head_A, worm_A : double
        Penalty for deviation from expected length.
    pobj, pneighbor, pbg : double, optional
        Probability of feature on object/neighbor/background.
    soft_min_size, hard_min_size = int, optional
        Minimum size of connected component to be considered as head initialization.
        If the currently best instantiation is above min_logpost then soft_min_size is used otherwise hard_min_size.
    min_logpost : double, optional
        See above.
    verbose : bool, optional
        Whether detection results should be printed and visualized.
        
    Returns
    -------
    coarse_points : list of pairs
        Detected worm instantiation (tail to head).
    log_posterior: double
        Log-posterior (relative to the background model) of the detection.
    head_length : int
        Length of the head instantiation.
    
    None is returned if no head features were found.
    """
    
    # set parameters for head detection
    _set_detection_parameters(head_features, head_lam, head_A, pobj, pneighbor, pbg, \
                              search_depth = 2, search_breadth = 30, remote = False)
    
    # find connected components of head features
    features = np.any(DETECTION_PARAMETERS.features, 2)    
    labels, nlabel = scipy.ndimage.label(features, structure=scipy.ndimage.generate_binary_structure(2,2))  # corresponds to a 9-element neighborhood
    if nlabel == 0: return None
    component_sizes = scipy.ndimage.sum(features, labels, range(nlabel+1))
    sorted_labels = np.argsort(component_sizes)[::-1]
    
    cdef:
        int *starting_orientations = [2, 6, 0, 4]  # down, up, right, left
        pair[int,int] starting_direction, zero_point, first_point, last_point
        Instantiation head_instantiation, best_head_instantiation
        Instantiation worm_instantiation_original, worm_instantiation_reverse
        int head_length, best_head_length, last_orientation
        int reverse_first_orientation, reverse_second_orientation, reverse_turn
        Instantiation best_worm_instantiation = Instantiation()
    
    # loop through connected components
    cdef int i, j, k
    for k in range(nlabel):
        
        # set parameters for head detection
        if k > 0:
            _set_detection_parameters(head_features, head_lam, head_A, pobj, pneighbor, pbg, \
                                      search_depth = 2, search_breadth = 30, remote = False)

        best_head_instantiation = Instantiation()

        # find corners
        starting_points = _find_border_points(features * (labels==sorted_labels[k]))  # top, bottom, left, right
    
        ### detect head ###
        
        for i in range(4):
            
            # create initial head instantiation
            head_instantiation = Instantiation()
            first_point = pair[int,int](starting_points[i][0], starting_points[i][1])
            head_instantiation.points.push_back(first_point)
            
            # block border points and neighbors in the reverse direction of the starting point
            blocked = np.zeros((DETECTION_PARAMETERS.features.shape[0],DETECTION_PARAMETERS.features.shape[1]), \
                               dtype=np.uint8)
            blocked[0,:] = True; blocked[DETECTION_PARAMETERS.features.shape[0]-1,:] = True
            blocked[:,0] = True; blocked[:,DETECTION_PARAMETERS.features.shape[1]-1] = True
            for j in range(3, 5+1): 
                blocked[first_point.first + DIRECTIONS[(starting_orientations[i]+j)%8].first, \
                        first_point.second + DIRECTIONS[(starting_orientations[i]+j)%8].second] = True
            #blocked[first_point.first, first_point.second] = True
            head_instantiation.blocked = blocked
            head_instantiation.last_orientation = starting_orientations[i]
            
            # compute initial likelihood and posterior
            starting_direction = DIRECTIONS[starting_orientations[i]]
            zero_point = pair[int,int](first_point.first-starting_direction.first, \
                                       first_point.second-starting_direction.second)
            head_instantiation.log_likelihood = _log_likelihood_change(zero_point, starting_orientations[i], blocked)
            head_instantiation.log_posterior = head_instantiation.log_likelihood - \
                                               DETECTION_PARAMETERS.A * (1 - DETECTION_PARAMETERS.lam)**2
    
            # store best posture
            head_instantiation, viz = _detect_posture1(head_instantiation, verbose)
            if head_instantiation.log_posterior > best_head_instantiation.log_posterior:
                best_head_instantiation = head_instantiation
                best_head_viz = viz
    
        # show result of head detection
        head_length = best_head_instantiation.points.size()
        
        ### detect body ###
        
        # set parameters for body detection
        _set_detection_parameters(worm_features, worm_lam, worm_A, pobj, pneighbor, pbg, \
                                  search_depth = 3, search_breadth = 100, remote = True)    
        
        # find body posture in original direction
        worm_instantiation_original = Instantiation(best_head_instantiation)
        worm_instantiation_original.log_likelihood = best_head_instantiation.log_posterior
        worm_instantiation_original.log_posterior = \
            worm_instantiation_original.log_likelihood - \
            DETECTION_PARAMETERS.A * (worm_instantiation_original.points.size() - DETECTION_PARAMETERS.lam)**2
        worm_instantiation_original, viz_orig = _detect_posture1(worm_instantiation_original, verbose)
            
        # adjust blocked points for reverse direction
        worm_instantiation_reverse = Instantiation(best_head_instantiation)
        last_point = best_head_instantiation.points.back()
        last_orientation = best_head_instantiation.last_orientation
        for i in range(-2,2+1):  # block neighbors of last point 
            worm_instantiation_reverse.blocked[last_point.first + DIRECTIONS[(last_orientation+i+8)%8].first, \
                                               last_point.second + DIRECTIONS[(last_orientation+i+8)%8].second] = True
        
        first_point = best_head_instantiation.points[0]
        reverse_first_orientation = ORIENTATIONS[pair[int,int](first_point.first - \
                                                                        best_head_instantiation.points[1].first, \
                                                                        first_point.second - \
                                                                        best_head_instantiation.points[1].second)]
        for i in range(-2,2+1):  # unblock neighbors of first point
            if 0 < first_point.first + DIRECTIONS[(reverse_first_orientation+i+8)%8].first \
                 < DETECTION_PARAMETERS.features.shape[0]-1 and \
               0 < first_point.second + DIRECTIONS[(reverse_first_orientation+i+8)%8].second \
               < DETECTION_PARAMETERS.features.shape[1]-1:   
                worm_instantiation_reverse.blocked[first_point.first + \
                                                   DIRECTIONS[(reverse_first_orientation+i+8)%8].first, \
                                                   first_point.second + \
                                                   DIRECTIONS[(reverse_first_orientation+i+8)%8].second] = False 
                                              
        # adjust points, last orientation and last turn for reverse direction
        worm_instantiation_reverse.points = best_head_instantiation.points[::-1]
        worm_instantiation_reverse.last_orientation = reverse_first_orientation
        reverse_second_orientation = ORIENTATIONS[pair[int,int](best_head_instantiation.points[1].first - \
                                                                         best_head_instantiation.points[2].first, \
                                                                         best_head_instantiation.points[1].second - \
                                                                         best_head_instantiation.points[2].second)]
        reverse_turn = (reverse_first_orientation - reverse_second_orientation + 8) % 8
        worm_instantiation_reverse.last_turn = -1 if reverse_turn == 7 else reverse_turn
    
        # find body posture in reverse direction
        worm_instantiation_reverse.log_likelihood = best_head_instantiation.log_posterior
        worm_instantiation_reverse.log_posterior = \
            worm_instantiation_reverse.log_likelihood - \
            DETECTION_PARAMETERS.A * (worm_instantiation_reverse.points.size() - DETECTION_PARAMETERS.lam)**2
        worm_instantiation_reverse, viz_reverse = _detect_posture1(worm_instantiation_reverse, verbose)
        
        # determine detection with higher posterior
        if worm_instantiation_original.log_posterior > best_worm_instantiation.log_posterior:
            best_worm_instantiation = worm_instantiation_original
            best_head_length = head_length
            if verbose:
                my_viz = np.empty((4,DETECTION_PARAMETERS.features.shape[0],DETECTION_PARAMETERS.features.shape[1],3), \
                                  dtype=np.uint8)
                my_viz[0,] = best_head_viz[3,]
                my_viz[1:4,] = viz_orig[1:4,]     
        if worm_instantiation_reverse.log_posterior > best_worm_instantiation.log_posterior:
            best_worm_instantiation = worm_instantiation_reverse
            best_head_length = head_length
            if verbose:
                my_viz = np.empty((4,DETECTION_PARAMETERS.features.shape[0],DETECTION_PARAMETERS.features.shape[1],3), \
                                  dtype=np.uint8)
                my_viz[0,] = best_head_viz[3,]
                my_viz[1:4,] = viz_reverse[1:4,]
    
        if (best_worm_instantiation.log_posterior > min_logpost \
            and component_sizes[sorted_labels[k+1]] < soft_min_size) \
        or component_sizes[sorted_labels[k+1]] < hard_min_size:
            break
        
    coarse_points = best_worm_instantiation.points
    
    if verbose:
        print("Found worm of length %d with head of length %d; log-lik = %.1f, log-posterior = %.1f \
               (initialized with %i clusters)" %(coarse_points.size(), best_head_length, \
                                                 best_worm_instantiation.log_likelihood, \
                                                 best_worm_instantiation.log_posterior, k+1))
        ##plot.show_image_collection(my_viz[:,1:-1,1:-1], ncol=4, rgb=True)
        
    return [(r-1,c-1) for (r,c) in coarse_points[::-1]], best_worm_instantiation.log_posterior, best_head_length  # inner feature array has one-point border


cdef _detect_posture2(Instantiation init_instantiation1, Instantiation init_instantiation2, bool doviz=False):
    """
    Runs a posture detection in two directions starting from the provided instantiations
    (and using the current detection setting).
    
    Parameters
    ----------
    init_instantiation1, init_instantiation2 : Instantiation
        Disjoint instantiations to start from.
    doviz : bool, optional
        Whether a visualization should be created.
        
    Returns
    -------
    best_instantiation : Instantiation
        The detected posture.
    viz : array
        Visualizations after 15, 30, 45 iterations, and the final detection.
    """

    # initialize best instantiations
    cdef Instantiation best_instantiation1 = init_instantiation1
    cdef Instantiation best_instantiation2 = init_instantiation2
    cdef int length = init_instantiation1.points.size() + init_instantiation2.points.size()
    cdef int last_improvement = length
    cdef double prior = - DETECTION_PARAMETERS.A * (length - DETECTION_PARAMETERS.lam)**2
    cdef double best_logpost = init_instantiation1.log_likelihood + init_instantiation2.log_likelihood + prior

    # initialize reference instantiations (similar to best_instantiation, but updated if likelihood increases, instead of posterior)
    cdef Instantiation reference_instantiation1  = init_instantiation1
    cdef Instantiation reference_instantiation2  = init_instantiation2

    # initialize candidate lists
    candidates1 = {init_instantiation1.points.back(): init_instantiation1}
    reference_candidates1 = candidates1
    candidates2 = {init_instantiation1.points.back(): init_instantiation2}
    reference_candidates2 = candidates2
    
    ### run detection (until no improvement happened for a while) ###
    
    cdef double loglik_change
    cdef double previous_length_best_loglik1 = init_instantiation1.log_likelihood
    cdef double previous_length_best_loglik2 = init_instantiation2.log_likelihood
    cdef Instantiation candidate, current_length_best_instantiation1, current_length_best_instantiation2
    cdef bool competition = False
    if doviz:
        visited = set(candidates1.keys()).union(set(candidates2.keys()))
        viz = np.zeros((4,DETECTION_PARAMETERS.features.shape[0],DETECTION_PARAMETERS.features.shape[1],3), \
                       dtype=np.uint8)
    else:
        viz = None
    
    while True:
        best_logpost
        # setup for next length
        length += 1
        new_candidates1 = {}
        current_length_best_instantiation1 = Instantiation()
        new_candidates2 = {}
        current_length_best_instantiation2 = Instantiation()

        # loop through candidates for first direction
        for candidate in candidates1.values():
            
            # if above the log-likelihood threshold...
            # ...or candidate is remote and above lower threshold (and candidate list not too long)
            if candidate.log_likelihood > previous_length_best_loglik1 - DETECTION_PARAMETERS.tau \
                or (DETECTION_PARAMETERS.remote \
                    and candidate.log_likelihood > previous_length_best_loglik1 - 2*DETECTION_PARAMETERS.tau \
                    and len(candidates1) <= DETECTION_PARAMETERS.search_breadth) \
                    and length - _intersection(best_instantiation1.points,candidate.points) \
                        >= DETECTION_PARAMETERS.lam//5:
                # make all moves and update the new best instantiation                
                current_length_best_instantiation1 = \
                    candidate._make_all_admissible_moves1(new_candidates1, current_length_best_instantiation1, \
                                                          reference_instantiation2.blocked)
        
        # loop through candidates for second direction
        for candidate in candidates2.values():
            
            # if above the log-likelihood threshold...
            # ...or candidate is remote and above lower threshold (and candidate list not too long)
            if candidate.log_likelihood > previous_length_best_loglik2 - DETECTION_PARAMETERS.tau \
                or (DETECTION_PARAMETERS.remote \
                    and candidate.log_likelihood > previous_length_best_loglik2 - 2*DETECTION_PARAMETERS.tau \
                    and len(candidates2) <= DETECTION_PARAMETERS.search_breadth) \
                    and length - _intersection(best_instantiation2.points,candidate.points) \
                        >= DETECTION_PARAMETERS.lam//5:
                # make all moves and update the new best instantiation                
                current_length_best_instantiation2 = \
                    candidate._make_all_admissible_moves1(new_candidates2, current_length_best_instantiation2, 
                                                          reference_instantiation1.blocked)
                
        # compute change in likelihoods (relative to the reference instantiations)
        loglik_change1 = current_length_best_instantiation1.log_likelihood - reference_instantiation1.log_likelihood
        loglik_change2 = current_length_best_instantiation2.log_likelihood - reference_instantiation2.log_likelihood
        prior = - DETECTION_PARAMETERS.A * (length - DETECTION_PARAMETERS.lam)**2

        if loglik_change1 >= loglik_change2:  # winner = 1
            
            # check if the best instantiation of this iteration has better posterior
            if reference_instantiation1.log_likelihood + loglik_change1 + reference_instantiation2.log_likelihood \
                                                                        + prior > best_logpost:
                best_logpost = reference_instantiation1.log_likelihood + loglik_change1 + \
                               reference_instantiation2.log_likelihood + prior
                best_instantiation1 = current_length_best_instantiation1
                best_instantiation2 = reference_instantiation2
                last_improvement = length
                assert length == best_instantiation1.points.size() + best_instantiation2.points.size()
            
            if loglik_change1 >= 0:  # end the competition (i.e., stop exploration)
                # update reference instantiation
                reference_instantiation1 = current_length_best_instantiation1
                reference_previous_length_best_loglik1 = current_length_best_instantiation1.log_likelihood
                reference_candidates1 = new_candidates1
                # reset looser
                previous_length_best_loglik2 = reference_previous_length_best_loglik2
                candidates2 = reference_candidates2
                competition = False
            
            else:  # start / continue competition (i.e., keep exploring)
                # backup candidates and current length best instantiation; 
                if not competition:
                    reference_candidates1 = candidates1
                    reference_previous_length_best_loglik1 = previous_length_best_loglik1
                    reference_candidates2 = candidates2
                    reference_previous_length_best_loglik2 = previous_length_best_loglik2
                    competition = True
                previous_length_best_loglik2 = current_length_best_instantiation2.log_likelihood
                candidates2 = new_candidates2
            
            # store results from current iteration    
            previous_length_best_loglik1 = current_length_best_instantiation1.log_likelihood
            candidates1 = new_candidates1

        else:  # winner = 2

            # check if the best instantiation of this iteration has better posterior
            if reference_instantiation1.log_likelihood + reference_instantiation2.log_likelihood + loglik_change2 \
                                                                        + prior > best_logpost:
                best_logpost = reference_instantiation1.log_likelihood + \
                               reference_instantiation2.log_likelihood + loglik_change2 + prior
                best_instantiation1 = reference_instantiation1
                best_instantiation2 = current_length_best_instantiation2
                last_improvement = length
                assert length == best_instantiation1.points.size() + best_instantiation2.points.size()
            
            if loglik_change2 >= 0:  # end the competition (i.e., stop exploration)
                # update reference instantiation
                reference_instantiation2 = current_length_best_instantiation2
                reference_previous_length_best_loglik2 = current_length_best_instantiation2.log_likelihood
                reference_candidates2 = new_candidates2
                # reset looser
                previous_length_best_loglik1 = reference_previous_length_best_loglik1
                candidates1 = reference_candidates1
                competition = False
            
            else:  # start / continue competition (i.e., keep exploring)
                # backup candidates and current length best instantiation; 
                if not competition:
                    reference_candidates1 = candidates1
                    reference_previous_length_best_loglik1 = previous_length_best_loglik1
                    reference_candidates2 = candidates2
                    reference_previous_length_best_loglik2 = previous_length_best_loglik2
                    competition = True
                previous_length_best_loglik1 = current_length_best_instantiation1.log_likelihood
                candidates1 = new_candidates1
            
            # store results from current iteration    
            previous_length_best_loglik2 = current_length_best_instantiation2.log_likelihood
            candidates2 = new_candidates2

        if doviz:
            visited = visited.union(set(new_candidates1.keys())).union(set(new_candidates2.keys()))
            viz_image = -1
            if length == 15:
                viz_image = 0
            elif length == 30:
                viz_image = 1
            elif length == 45:
                viz_image = 2
            if viz_image != -1:
                rr, cc = zip(*visited)
                viz[viz_image, rr, cc] = (0,128,128)
                rr, cc = np.where(np.any(DETECTION_PARAMETERS.features,2))
                viz[viz_image, rr,cc,] = 128
                rr, cc = zip(*candidates1.keys())
                viz[viz_image, rr, cc] = (128,0,128)
                rr, cc = zip(*candidates2.keys())
                viz[viz_image, rr, cc] = (128,0,128)
                rr, cc = zip(*best_instantiation1.points)
                viz[viz_image, rr,cc,] = 255
                rr, cc = zip(*best_instantiation2.points)
                viz[viz_image, rr,cc,] = 255
        
        # stop search if no improvement has happened for too long
        if length == last_improvement + DETECTION_PARAMETERS.search_depth:
            break
    
    # create full instantiation
    cdef Instantiation worm_instantiation = Instantiation()
    worm_instantiation.log_likelihood = best_instantiation1.log_likelihood + best_instantiation2.log_likelihood
    worm_instantiation.log_posterior = best_logpost
    worm_instantiation.points = best_instantiation1.points[::-1] + best_instantiation2.points[:]
    if doviz:
        rr, cc = zip(*visited)
        viz[3, rr, cc] = (0,128,128)
        rr, cc = np.where(np.any(DETECTION_PARAMETERS.features,2))
        viz[3, rr,cc,] = 128
        rr, cc = zip(*worm_instantiation.points)
        viz[3, rr,cc,] = 255
    
    return worm_instantiation, viz


def detect_worm2(pair[int,int] starting_point, int starting_orientation, \
                 cnp.ndarray[cnp.uint8_t, ndim=3, mode='c'] head_features, double head_lam, double head_A, \
                 cnp.ndarray[cnp.uint8_t, ndim=3, mode='c'] worm_features, double worm_lam, double worm_A, \
                 double pobj = .7, double pneighbor = .5, double pbg = .3, bool verbose=False):
    """
    Detects the worm posture by starting a search at the provided point going in both directions. 
    
    
    Parameters
    ----------
    starting_point : pair[int,int]
        Starting point of the search.
    starting_orientation : int
        Initial search orientation. 
    head_features, worm_features : array
        Detected head/worm features.
    head_lam, worm_lam : double
        Expected length of head/worm posture (in coarse blocks). 
    head_A, worm_A : double
        Penalty for deviation from expected length.
    pobj, pneighbor, pbg : double, optional
        Probability of feature on object/neighbor/background.
    verbose : bool, optional
        Whether detection results should be printed and visualized.
        
    Returns
    -------
    coarse_points : list of pairs
        Detected worm instantiation (tail to head).
    log_posterior: double
        Log-posterior (relative to the background model) of the detection.
    head_length : int
        Length of the head instantiation.
    """
    
    starting_point.first += 1
    starting_point.second += 1
    
    cdef pair[int,int] starting_direction = DIRECTIONS[starting_orientation]
    cdef pair[int,int] starting_point1 = pair[int,int](starting_point.first + starting_direction.first, \
                                                       starting_point.second + starting_direction.second)
    cdef pair[int,int] starting_point2 = starting_point
    cdef int reverse_starting_orientation = (starting_orientation + 4) % 8
    
    # set parameters for worm detection
    _set_detection_parameters(worm_features, worm_lam, worm_A, pobj, pneighbor, pbg, \
                              search_depth = 3, search_breadth = 100, remote = True) 
    
    # initialize instantiations
    cdef Instantiation init_instantiation1 = Instantiation()
    init_instantiation1.points.push_back(starting_point1)
    init_instantiation1.last_orientation = starting_orientation
    cdef Instantiation init_instantiation2 = Instantiation()
    init_instantiation2.points.push_back(starting_point2)
    init_instantiation2.last_orientation = reverse_starting_orientation
    
    # block border points and neighbors in the reverse direction of the starting point
    blocked1 = np.zeros((DETECTION_PARAMETERS.features.shape[0],DETECTION_PARAMETERS.features.shape[1]), dtype=np.uint8)
    blocked1[0,:] = True; blocked1[DETECTION_PARAMETERS.features.shape[0]-1,:] = True
    blocked2 = np.zeros((DETECTION_PARAMETERS.features.shape[0],DETECTION_PARAMETERS.features.shape[1]), dtype=np.uint8)
    blocked2[0,:] = True; blocked2[DETECTION_PARAMETERS.features.shape[0]-1,:] = True
    blocked1[:,0] = True; blocked1[:,DETECTION_PARAMETERS.features.shape[1]-1] = True
    for j in range(3, 5+1): 
        blocked1[starting_point1.first + DIRECTIONS[(starting_orientation+j)%8].first, \
                 starting_point1.second + DIRECTIONS[(starting_orientation+j)%8].second] = True
    ##blocked1[starting_point1.first, starting_point1.second] = True
    init_instantiation1.blocked = blocked1
    blocked2[:,0] = True; blocked2[:,DETECTION_PARAMETERS.features.shape[1]-1] = True
    for j in range(3, 5+1): 
        blocked2[starting_point2.first + DIRECTIONS[(reverse_starting_orientation+j)%8].first, \
                 starting_point2.second + DIRECTIONS[(reverse_starting_orientation+j)%8].second] = True
    ##blocked2[starting_point2.first, starting_point2.second] = True
    init_instantiation2.blocked = blocked2
    
    # compute initial likelihood and posterior
    init_instantiation1.log_likelihood = _log_likelihood_change(starting_point2, starting_orientation, blocked1)
    init_instantiation1.log_posterior = init_instantiation1.log_likelihood - \
                                        DETECTION_PARAMETERS.A * (1 - DETECTION_PARAMETERS.lam)**2
    init_instantiation2.log_likelihood = _log_likelihood_change(starting_point1, reverse_starting_orientation, blocked2)
    init_instantiation2.log_posterior = init_instantiation2.log_likelihood - \
                                        DETECTION_PARAMETERS.A * (1 - DETECTION_PARAMETERS.lam)**2

    # find worm posture
    cdef Instantiation worm_instantiation
    worm_instantiation, viz = _detect_posture2(init_instantiation1, init_instantiation2, verbose)
    
    ### determine head ending ###
    
    # set parameters for head detection
    _set_detection_parameters(head_features, head_lam, head_A, pobj, pneighbor, pbg, \
                              search_depth = 3, search_breadth = 100, remote = True)
    
    cdef int worm_length = worm_instantiation.points.size()
    not_blocked = np.zeros((DETECTION_PARAMETERS.features.shape[0],DETECTION_PARAMETERS.features.shape[1]), \
                           dtype=np.uint8)
    
    # find best head instantiation from one ending
    cdef int orientation = ORIENTATIONS[pair[int,int](worm_instantiation.points[1].first - \
                                                      worm_instantiation.points[0].first, \
                                                      worm_instantiation.points[1].second - \
                                                      worm_instantiation.points[0].second)]
    cdef pair[int,int] last_point = pair[int,int](worm_instantiation.points[0].first - DIRECTIONS[orientation].first, \
                                                  worm_instantiation.points[0].second - DIRECTIONS[orientation].second)
    cdef int length = 1, head_length
    cdef double loglik = 0, logpost, best_loglik = -INFINITY, best_logpost = -INFINITY
    while True:
        # compute head posterior
        loglik += _log_likelihood_change(last_point, orientation, not_blocked)
        logpost = loglik - DETECTION_PARAMETERS.A * (length - DETECTION_PARAMETERS.lam)**2
        if logpost > best_logpost:
            head_length = length
            best_loglik = loglik
            best_logpost = logpost
        
        # update last point and orientation
        length += 1
        if length > 2*head_lam or length > worm_length:
            break
        orientation = ORIENTATIONS[pair[int,int](worm_instantiation.points[length-1].first - \
                                                 worm_instantiation.points[length-2].first, \
                                                 worm_instantiation.points[length-1].second - \
                                                 worm_instantiation.points[length-2].second)]
        last_point = worm_instantiation.points[length-2]

    # find best head instantiation from the other ending
    orientation = ORIENTATIONS[pair[int,int](worm_instantiation.points[worm_length-2].first - \
                                             worm_instantiation.points[worm_length-1].first, \
                                             worm_instantiation.points[worm_length-2].second - \
                                             worm_instantiation.points[worm_length-1].second)]
    last_point = pair[int,int](worm_instantiation.points[worm_length-1].first - DIRECTIONS[orientation].first, \
                               worm_instantiation.points[worm_length-1].second - DIRECTIONS[orientation].second)
    length = 1
    loglik = 0
    while True:
        # compute head posterior
        loglik += _log_likelihood_change(last_point, orientation, not_blocked)
        logpost = loglik - DETECTION_PARAMETERS.A * (length - DETECTION_PARAMETERS.lam)**2
        if logpost > best_logpost:
            head_length = -length
            best_loglik = loglik
            best_logpost = logpost
        
        # update last point and orientation
        length += 1
        if length > 2*head_lam or length > worm_length:
            break
        orientation = ORIENTATIONS[pair[int,int](worm_instantiation.points[worm_length-length].first - \
                                                 worm_instantiation.points[worm_length-length+1].first, \
                                                 worm_instantiation.points[worm_length-length].second - \
                                                 worm_instantiation.points[worm_length-length+1].second)]
        last_point = worm_instantiation.points[worm_length-length+1]
    
    if head_length < 0:
        worm_instantiation.points = worm_instantiation.points[::-1]
        head_length *= -1

    ### determine total likelihood ###
    
    # set parameters for worm detection
    _set_detection_parameters(worm_features, worm_lam, worm_A, pobj, pneighbor, pbg, \
                              search_depth = 3, search_breadth = 100, remote = True)

    # compute body likelihood
    loglik = best_loglik
    for length in range(head_length, worm_length):
        last_point = worm_instantiation.points[length]
        orientation = ORIENTATIONS[pair[int,int](worm_instantiation.points[length+1].first - last_point.first, \
                                                 worm_instantiation.points[length+1].second - last_point.second)]
        loglik += _log_likelihood_change(last_point, orientation, not_blocked)
    logpost = loglik - DETECTION_PARAMETERS.A * (worm_length - DETECTION_PARAMETERS.lam)**2
        
    # show result of worm detection
    if verbose:
        print("Found worm of length %d with head of length %d; log-lik = %.1f, log-posterior = %.1f" \
              %(worm_instantiation.points.size(), head_length, loglik, logpost))
        ##plot.show_image_collection(viz[:,1:-1,1:-1], rgb=True)
    
    return [(r-1,c-1) for (r,c) in worm_instantiation.points[::-1]], logpost, head_length  # inner feature array has one-point border
