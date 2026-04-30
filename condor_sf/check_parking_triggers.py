import uproot
path = 'root://cmsxrootd.fnal.gov///store/data/Run2023C/ParkingVBF0/NANOAOD/22Sep2023_v3-v1/2550000/27f82e00-fdb5-4f73-a8ba-637274804d27.root'
f = uproot.open(path)
t = f['Events']
branches = t.keys()
muon_trigs = sorted([b for b in branches if b.startswith('HLT_') and 'Mu' in b])
vbf_trigs  = sorted([b for b in branches if b.startswith('HLT_') and 'VBF' in b])
print('=== MUON TRIGGERS ===')
print('\n'.join(muon_trigs))
print('\n=== VBF TRIGGERS ===')
print('\n'.join(vbf_trigs))
