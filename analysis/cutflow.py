
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import matplotlib.pyplot as plt
from itertools import combinations
from scipy import constants
from math import acos
import visualization
import physics
import ROOT as root
import numpy as np
import detector
import util
import time
import scipy.optimize as spo
import sys
import glob

ncuts = 7



class sample_space():


    def __init__(self,file):

        self.file = file


    def inside_box(self, x, y, z):

        # box_lims = [[-5000., 5000.], [5900.0, 9000.0], [6900., 17100.]]
        box_lims = [[-5000., 5000.], [6000.0, 8917.0], [7000., 17000.]]

        if box_lims[0][0] < x and box_lims[0][1] > x:
            if box_lims[1][0] < y and box_lims[1][1] > y:
                if box_lims[2][0] < z and box_lims[2][1] > z:
                    return True

        return False


    def in_layer(self,y_val):

        #LAYERS_Y=[[6001.0, 6004.0],  [6104.0, 6107.0], [8001.0, 8004.0], [8104.0, 8107.0], \
        #[8501.0, 8504.0], [8604.0, 8607.0], [8707.0, 8710.0], [8810.0, 8813.0], [8913.0, 8916.0]]

        LAYERS_Y=[[6003.0 + 547, 6006.0 + 547],  #layer 0 (floor)
		[6106.0 + 547, 6109.0 + 547], #layer 1 (floor)
		[6209.0 + 547, 6212.0 + 547], #layer 2 (floor)
		[8003.0 + 547, 8006.0 + 547], #layer 3
		[8106.0 + 547, 8109.0 + 547], #layer 4
		[8503.0 + 547, 8506.0 + 547], #layer 5
		[8606.0 + 547, 8609.0 + 547], #layer 6
		[8709.0 + 547, 8712.0 + 547], #layer 7
		[8812.0 + 547, 8815.0 + 547], #layer 8
		[8915.0 + 547, 8918.0 + 547]] #layer 9

        for n in range(len(LAYERS_Y)):
            _min = LAYERS_Y[n][0]
            _max = LAYERS_Y[n][1]
            if (y_val > _min) and (y_val < _max):
                    return n

        return 999


    def get_states_m(self):

        global ncuts

        tracking_file = root.TFile.Open(self.file)
        self.tree = tracking_file.Get("integral_tree")
        self.tree.SetBranchStatus("*", 0)

        self.tree.SetBranchStatus("Digi_y", 1)
        self.tree.SetBranchStatus("Digi_time", 1)
        self.tree.SetBranchStatus("NumTracks_k_m", 1)

        self.tree.SetBranchStatus("Vertex_k_m_t", 1)
        self.tree.SetBranchStatus("Vertex_k_m_x", 1)
        self.tree.SetBranchStatus("Vertex_k_m_y", 1)
        self.tree.SetBranchStatus("Vertex_k_m_z", 1)
        self.tree.SetBranchStatus("Vertex_k_m_trackIndices", 1)
        self.tree.SetBranchStatus("NumVertices_k_m", 1)

        self.tree.SetBranchStatus("Track_k_m_x0", 1)
        self.tree.SetBranchStatus("Track_k_m_velX", 1)
        self.tree.SetBranchStatus("Track_k_m_velY", 1)
        self.tree.SetBranchStatus("Track_k_m_velZ", 1)
        self.tree.SetBranchStatus("Track_k_m_hitIndices", 1)
        self.tree.SetBranchStatus("Track_k_m_expected_hit_layer", 1)

        cut_vectors = []

        for event_number in range(int(self.tree.GetEntries())):

            current_vector = np.zeros(1+ncuts) # first will be the event number

            self.tree.GetEntry(event_number)
            if event_number % 500 == 0:
                print("event:", event_number)

            # record the event number for reference
            current_vector[0] = event_number

            # cut on the number of made tracks in the event
            current_vector[1] = self.tree.NumTracks_k_m

            # cut on the number of vertices made in the event
            current_vector[2] = len(self.tree.Vertex_k_m_x)

            inside = False

            for num in range(len(self.tree.Vertex_k_m_x)):
                vtxx, vtxy, vtxz = self.tree.Vertex_k_m_x[num], self.tree.Vertex_k_m_y[num], self.tree.Vertex_k_m_z[num]
                if self.inside_box(vtxx, vtxy, vtxz):
                    inside = True

            # check if at least one vertex is in the detector
            current_vector[3] = int(inside) # fiducial

            floorveto = False
            for num in range(len(self.tree.Vertex_k_m_x)): # loop over vertices
                vertexveto = False
                for hit in range(len(self.tree.Digi_y)): # loop over the hits
                    if self.in_layer(self.tree.Digi_y[hit]) <= 2: # check if layer index is a floor layer index
                        if self.tree.Vertex_k_m_t[num] > self.tree.Digi_time[hit]: # check if the floor hit happened before the vertex
                            vertexveto = True
                if not vertexveto: # There's a vertex in the event that isn't vetoed
                    break
            if len(self.tree.Vertex_k_m_x) != 0:
                if vertexveto:
                    floorveto = True

            # track floor veto
            current_vector[4] = int(not floorveto) # cut if floorveto == True => 0 == current_vector[4] (< 1 is true)


            expectedveto = False
            fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)
            fullhitindices = util.unzip(self.tree.Track_k_m_hitIndices)
            fullexp_layers = util.unzip(self.tree.Track_k_m_expected_hit_layer)

            for j in range(int(self.tree.NumVertices_k_m)):
                bottomlayer_exp = []
                bottomlayer_hits = []
                trackindices = fulltrackindices[j]

                for track in trackindices:
                    exp_layers = fullexp_layers[int(track)]  # layers the track is expected to hit

                    for hit in exp_layers:
                        if hit < 2:
                            bottomlayer_exp.append(hit) # a floor layer hit is expected, record the layer

                    hitindices = fullhitindices[int(track)]

                    for i in hitindices:
                        hity = self.tree.Digi_y[i]

                        if self.in_layer(hity) <= 2:
                                bottomlayer_hits.append(hity) # hit is in a floor layer

                # we get 3 or more expected hits, but no hits, veto
                if (len(bottomlayer_hits) < 1 and len(bottomlayer_exp) >= 3):
                    expectedveto = True

            current_vector[5] = int(expectedveto) # bottom layer hits

            open_angles = []
            fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)
            for k1 in range(len(self.tree.Vertex_k_m_x)):
                trackindices = fulltrackindices[k1]
                combolist = list(combinations(trackindices, 2))

                for combo in combolist:
                    combo = [int(combo[0]),int(combo[1])]

                    cos_opening_angle = self.tree.Track_k_m_velX[combo[0]]*self.tree.Track_k_m_velX[combo[1]] + \
                                        self.tree.Track_k_m_velY[combo[0]]*self.tree.Track_k_m_velY[combo[1]] + \
                                        self.tree.Track_k_m_velZ[combo[0]]*self.tree.Track_k_m_velZ[combo[1]]
                    mag1 = np.sqrt(self.tree.Track_k_m_velX[combo[0]]**2 + self.tree.Track_k_m_velY[combo[0]]**2 + self.tree.Track_k_m_velZ[combo[0]]**2)
                    mag2 = np.sqrt(self.tree.Track_k_m_velX[combo[1]]**2 + self.tree.Track_k_m_velY[combo[1]]**2 + self.tree.Track_k_m_velZ[combo[1]]**2)
                    cos_opening_angle = cos_opening_angle/( mag1*mag2 )
                    open_angles.append(acos(cos_opening_angle))

            # sort in terms of descending opening angle
            open_angles.sort(reverse = True)

            if len(open_angles) != 0:
                        # cut on largest opening angle
                	current_vector[6] = open_angles[0]


            # check if any vertices are below all corresponding track hits
            topoveto = 0

            fulltrackindices = util.unzip(self.tree.Vertex_k_m_trackIndices)
            fullhitindices = util.unzip(self.tree.Track_k_m_hitIndices)

            for k1 in range(len(self.tree.Vertex_k_m_x)):
                vertexveto = False
                trackindices = fulltrackindices[k1]
                for t1 in range(len(self.tree.Track_k_m_x0)):
                    hitindices = fullhitindices[t1]
                    for hit in hitindices:
                        if self.tree.Digi_y[hit] < self.tree.Vertex_k_m_y[k1]:
                            vertexveto = True
                if vertexveto:
                    topoveto += 1

            if topoveto < current_vector[2]: # vetoed vertices < total vertices
                # see if there is at least one vertex below all hits in its made tracks
                current_vector[7] = 1

            cut_vectors.append(current_vector)

        self.cut_vectors = np.array(cut_vectors)

        return self.cut_vectors




class scissors():
    '''since this is what does the cutting!'''

    def __init__(self,file):

        self.file = file

    def store_space(self):

        self.space = sample_space(self.file)

        # self.cut_vectors = self.space.get_states()
        self.cut_vectors = self.space.get_states_m()

        print('number of events is: ',len(self.cut_vectors))


    def cut(self, p0, p1, p2, p3, p4, p5, p6):

        global ncuts

        self.cuts = [p0, p1, p2, p3, p4, p5, p6]

        self.flows = np.zeros(ncuts)
        local_vectors = np.copy(self.cut_vectors)
        self.survivor_inds = []

        for i in range(ncuts):

            inds = np.where(local_vectors[:,i+1] < self.cuts[i]) # +1 since first is ev number
            local_vectors = np.delete(local_vectors, inds, 0)

            self.flows[i] += np.shape(local_vectors)[0]
            self.survivor_inds.append(local_vectors[:,0])

        self.survivors = self.flows[-1]

#        print("flows are ",self.flows)
#        print("with parameters ",self.cuts)



class par_func():
    ''' handle multiple files to cut and create a function that takes cut parameters and returns signal to background ratios'''

    def __init__(self,backgrounds,signals):

        if type(backgrounds) == list and type(signals) == list:
            self.b_files = backgrounds
            self.s_files = signals

        else:
            self.b_files = [backgrounds]
            self.s_files = [signals]

    def get_states(self):

        self.drawers = [[],[]]

        i = 0

        for files in [self.s_files, self.b_files]:
            for fl in files:

                print(fl)

                scissor = scissors(fl)

                scissor.store_space()
                self.drawers[i].append(scissor)

            i += 1


#    def cut(self, p0, p1, p2, p3, p4, p5):
    def cut(self, p0, p1, p4):
        ''' calculate and return signal to background ratio for all handled files'''

        i = 0
        self.scores = [[],[]]

        for drawer in self.drawers:
            for scissor in drawer:

#                scissor.cut(p0, p1, p2, p3, p4, p5)
                scissor.cut(p0, p1, 1, 1, 1, p4, 1)
                self.scores[i].append(scissor.survivors / \
				      np.shape(scissor.cut_vectors)[0])

            i += 1

        return np.sum(self.scores[0]) if np.sum(self.scores[1]) == 0 \
                                      else np.sum(self.scores[0]) / np.sum(self.scores[1])


#    def get_flows(self, p0, p1, p2, p3, p4, p5):
    def get_flows(self, p0, p1, p4):

       flows = [[],[]]

       i = 0

       for drawer in self.drawers:
            for scissor in drawer:

#                scissor.cut(p0, p1, p2, p3, p4, p5)
                scissor.cut(p0, p1, 1, 1, 1, p4, 1)
                flows[i].append(scissor.flows)

            i += 1

       return flows



def main():

    directory = sys.argv[1]

    b_files = [filename for filename in glob.iglob(directory+'stat_0_*.root', recursive=True)]
    s_files = [filename for filename in glob.iglob(directory+'stat_2_*.root', recursive=True)]

    print(b_files)
    print(s_files)

    for i in range(len(b_files)):

        func = par_func(b_files[i],s_files[i])

        func.get_states()

        '''
#       bounds = [(-0.01,5.01),(-0.01,1.01),(-0.01,1.01),(-0.01,1.01),(-1.01,1.01),(-0.01,1.01)]
        bounds = [(-0.01,5.01),(-0.01,3.01),(-1.01,1.01)]

        opts = ({"minimize_every_iter":True,
                 "disp":False})

        # need the lambda function since shgo doesn't like methods!
        x = spo.shgo(lambda X: -func.cut(*X), bounds, iters=6, options=opts)

        print("pars are ",x)
        print("optimal is ",func.cut(*x['x']))

        flows = func.get_flows(*x['x'])
        print("signal flows are \n ",flows[0],"\nbackground flows are \n",flows[1])
        '''
        flows = func.get_flows(2,1,-1)
        print("optimal is ",func.cut(2,1,-1))
        print("signal flows are \n ",flows[0],"\nbackground flows are \n",flows[1])




if __name__ == '__main__':

    '''
    main()

    '''
    scissor = scissors(sys.argv[1])

    scissor.store_space()
    scissor.cut(2,1,1,1,-1,-1,1)

    print("flows are ",scissor.flows)
#    print("events with vertices ",scissor.survivor_inds[1])
#    print("survivor inds are ",np.array(scissor.survivor_inds[6]))