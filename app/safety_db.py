BANNED_PESTICIDES = {
    "2_4_d": {
        "name": "2,4-D",
        "smiles": "CC(=O)OC1=CC(=C(C=C1Cl)O)Cl",
    },
    "acifluorfen_sodium": {
        "name": "Acifluorfen (neutral form)",
        # sodium salt often written with [Na+], but neutral is fine for similarity
        "smiles": "O=C(O)C1=CC(=C(C=C1[N+](=O)[O-])OCC2=CC=CC=C2Cl)C(F)(F)F",
    },
    "biphenyl": {
        "name": "Biphenyl",
        "smiles": "C1=CC=C(C=C1)C2=CC=CC=C2",
    },
    "bromophos_ethyl": {
        "name": "Bromophos-ethyl",
        "smiles": "CCOP(=O)(OCC)OC1=NC(=C(C(=C1)Cl)Cl)Br",
    },
    "captan": {
        "name": "Captan",
        "smiles": "CC1(C(=O)N(C(=O)N1C(Cl)Cl)C(=O)N(CCl)CCl)C(=O)N(CCl)CCl",
    },
    "chlorpropham": {
        "name": "Chlorpropham",
        "smiles": "CCOC(=O)NC1=CC(=C(C=C1)Cl)C",
    },
    "cholecalciferol": {
        "name": "Cholecalciferol (Vitamin D3)",
        "smiles": "CC(C)CCCC(C)CCCC(C)C1CCC2C1(CCCC2=CC=C3C(=O)CCC(=C3)C)C",
    },
    "cyanazine": {
        "name": "Cyanazine",
        "smiles": "CCN(C)C1=NC(=NC(=N1)Cl)NC(C)C",
    },
    "demeton_methyl": {
        "name": "Demeton-methyl",
        "smiles": "COCSP(=O)(OC)OCCSC",
    },
    "dichlobenil": {
        "name": "Dichlobenil",
        "smiles": "C1=CC(=CC(=C1Cl)Cl)C#N",
    },
    "dichlorprop": {
        "name": "Dichlorprop",
        "smiles": "CC(C(=O)O)C1=CC(=CC=C1O)Cl",
    },
    "dodine": {
        "name": "Dodine",
        "smiles": "CCCCN(CCCN(CCCC)CCCC)C(=O)NC1=CC=CC=C1",
    },
    "ethion": {
        "name": "Ethion",
        "smiles": "CCOP(=S)(OCC)SCCSP(=S)(OCC)OCC",
    },
    "fenbuconazole": {
        "name": "Fenbuconazole",
        "smiles": "CC1=CC(=NO1)C2=CN=CC(=C2Cl)C3=CC=C(C=C3)Cl",
    },
    "fenhexamid": {
        "name": "Fenhexamid",
        "smiles": "CC(C)(C1=CC(=CC(=C1O)Cl)C(=O)NCCOC2=CC=CC=C2)O",
    },
    "ferbam": {
        "name": "Ferbam",
        "smiles": "CNC(=S)N(CCNC(=S)N)C(=S)N",
    },
    "fluazinam": {
        "name": "Fluazinam",
        "smiles": "CC1=CC(=NN1C2=NC(=CC(=C2Cl)Cl)C(F)(F)F)C(F)(F)F",
    },
    "flusulfamide": {
        "name": "Flusulfamide",
        "smiles": "CC1=CC(=NN1C2=NC(=CC(=C2F)S(=O)(=O)N)F)C(F)(F)F",
    },
    "fluvalinate": {
        "name": "Fluvalinate",
        "smiles": "CC(C)C(C(=O)OC(C#N)C1=CC(=CC=C1)OC2=CC=CC=C2)NC3=CC=C(C=C3)C(F)(F)F",
    },
    "forchlorfenuron": {
        "name": "Forchlorfenuron",
        "smiles": "CC1=NC(=NC(=C1)NC(=O)N)C2=CC=CC=C2Cl",
    },
    "furfural": {
        "name": "Furfural",
        "smiles": "O=CC1=CC=CO1",
    },
    "halosulfuron_methyl": {
        "name": "Halosulfuron-methyl",
        "smiles": "COC(=O)N1C(=O)N(S(=O)(=O)C2=CC(=CC=C21)Cl)C3=NC(=NC=C3)N",
    },
    "imazalil": {
        "name": "Imazalil (base)",
        "smiles": "CC1=CC(=NO1)C2=CC=CC(=C2Cl)CCN3CCN(CC3)CC",
    },
    "lactofen": {
        "name": "Lactofen",
        "smiles": "COC1=CC=C(C=C1OC(=O)C2=CC(=C(C=C2)Cl)[N+](=O)[O-])C(C)C",
    },
    "mecoprop": {
        "name": "Mecoprop (MCPP)",
        "smiles": "CC(C(=O)O)C1=CC(=CC=C1O)Cl",
    },
    "meptyldinocap": {
        "name": "Meptyldinocap",
        "smiles": "CCOC(=O)C1=CC(=CC(=C1OCC(C)C)C(C)C)C(C)C",
    },
    "pyrethrins": {
        "name": "Pyrethrins (represented as Pyrethrin I)",
        "smiles": "CC1=C(C(=O)C[C@@H]1OC(=O)[C@@H]2[C@H](C2(C)C)C=C(C)C)C/C=C\\C=C",
    },
    "pyrimidifen": {
        "name": "Pyrimidifen",
        "smiles": "CC1=NC(=NC(=C1SC2=NC(=NC(=C2)C3=CC=CC=C3)C)C)C4=CC=CC=C4",
    },
    "simazine": {
        "name": "Simazine",
        "smiles": "CCNC1=NC(=NC(=N1)NCC)Cl",
    },
    "tau_fluvalinate": {
        "name": "Tau-fluvalinate",
        "smiles": "CC(C)[C@@H](Nc1ccc(cc1Cl)C(F)(F)F)C(=O)OC(C#N)c2cccc(OC3=CC=CC=C3)c2",
    },
    "tebuconazole": {
        "name": "Tebuconazole",
        "smiles": "CC(C)(C)C(O)(CCC1=CC=C(Cl)C=C1)CN2C=NC=N2",
    },
    "thiabendazole": {
        "name": "Thiabendazole",
        "smiles": "c1ccc2[nH]c(nc2c1)-c3cscn3",
    },
    "tribufos": {
        "name": "Tribufos",
        "smiles": "O=P(SCCCC)(SCCCC)SCCCC",
    },
    "trichloroacetic_acid": {
        "name": "Trichloroacetic acid",
        "smiles": "O=C(O)C(Cl)(Cl)Cl",
    },
}
