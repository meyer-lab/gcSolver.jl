



# Assert the conservation of species throughout the experiment.
function assertConservation(y, y0, IDX)
        [1, 4, 5, 7, 8, 11, 12, 14, 15]),  # IL2Rb
    np.array([0, 3, 5, 6, 8]),  # IL2Ra
    np.array([9, 10, 12, 13, 15]),  # IL15Ra
    np.array([16, 17, 18]),  # IL7Ra
    np.array([19, 20, 21]),  # IL9R
    np.array([22, 23, 24]),  # IL4Ra
    np.array([25, 26, 27]),  # IL21Ra
    np.array([2, 6, 7, 8, 13, 14, 15, 18, 21, 24, 27]), #gc

        species_delta = y - y0

        # Check for conservation of species sum
        self.assertAlmostEqual(np.sum(species_delta[IDX]), 0.0, msg=str(IDX))
end



###
# In the absence of trafficking, mass balance should hold in both compartments.
rxntfR = self.rxntfR.copy()
rxntfR[17:30] = 0.0

dy = gcSolver.fullModel(y0, 0.0, rxntfR)

# Check for conservation of each surface receptor
assertConservation(dy, 0.0, idxs)

# Check for conservation of each endosomal receptor
assertConservation(dy, 0.0, idxs + 28)


