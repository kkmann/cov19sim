dm = HeavyTailsModel(LarremoreModel(0.02), l = 3.0, scale = 1.0, df = 3.0)

vl = log10.(sample(dm).vl)
