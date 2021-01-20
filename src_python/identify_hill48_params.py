import csv

with open("Datas/AA2090-T3.csv") as f:
    reader = csv.DictReader(f)
    exp_data = [row for row in reader]


for data in exp_data:
    orientation = data["orientation"]
    yield_stress = float(data["normalized_yield_stress"])
    r_value = float(data["r_value"])
    if orientation == "0":
        s0 = yield_stress
        r0 = r_value
    elif orientation == "45":
        s45 = yield_stress
        r45 = r_value
    elif orientation == "90":
        s90 = yield_stress
        r90 = r_value
    elif orientation == "EB":
        sb = yield_stress
        rb = r_value


stressF = (1/s90**2 + 1/sb**2 - 1/s0**2)/2
stressG = (1/sb**2 + 1/s0**2 - 1/s90**2)/2
stressH = (1/s0**2 + 1/s90**2 - 1/sb**2)/2
stressN = ((2/s45)**2 - (1/sb)**2)/2

print("stressF", stressF)
print("stressG", stressG)
print("stressH", stressH)
print("stressN", stressN)

rF = r0/(r90*(1.0+r0))
rG = 1.0/(1.0+r0)
rH = r0/(1.0+r0)
rN = 0.5*(r0+r90)*(1.0+2.0*r45)/(r90*(1.0+r0))
print("rF", rF)
print("rG", rG)
print("rH", rH)
print("rN", rN)
