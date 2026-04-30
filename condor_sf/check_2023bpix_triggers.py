import uproot

path = 'root://cmsxrootd.fnal.gov//store/data/Run2023D/Muon0/NANOAOD/22Sep2023_v1-v1/2530000/dd9c19c9-84b4-47e1-bc42-34320f55faba.root'
f = uproot.open(path)
t = f['Events']
branches = t.keys()

# What we currently use for 2023BPix in the code
hbb_triggers_2023bpix = [
    'AK8PFJet250_SoftDropMass40_PNetBB0p06',
    'PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70',
    'VBF_DiPFJet125_45_Mjj720_Detajj3p0'
]

print('=== CHECKING CONFIGURED HBB TRIGGERS ===')
for trig in hbb_triggers_2023bpix:
    present = ('HLT_' + trig) in branches
    print(f"  HLT_{trig}: {'FOUND' if present else 'MISSING'}")

print('\n=== ALL AK8 JET TRIGGERS ===')
ak8_trigs = sorted([b for b in branches if b.startswith('HLT_') and 'AK8' in b])
for t in ak8_trigs:
    print(f'  {t}')

print('\n=== ALL PFHT QUAD JET TRIGGERS ===')
quad_trigs = sorted([b for b in branches if b.startswith('HLT_') and 'QuadPFJet' in b])
for t in quad_trigs:
    print(f'  {t}')

print('\n=== ALL PNet/ParticleNet TRIGGERS ===')
pnet_trigs = sorted([b for b in branches if b.startswith('HLT_') and ('PNet' in b or 'ParticleNet' in b)])
for t in pnet_trigs:
    print(f'  {t}')

print('\n=== ALL PFHT TRIGGERS ===')
pfht_trigs = sorted([b for b in branches if b.startswith('HLT_') and 'PFHT' in b])
for t in pfht_trigs:
    print(f'  {t}')
