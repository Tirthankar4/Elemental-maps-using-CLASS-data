############ For highland detection ############
"""
Define the function to return true if the given point is inside the highland region
"""
def highland1(lat, lon):
    # Highland 1
    highlandName = "Highland 1"
    if lat >= -20 and lat <= 20 and lon >= -40 and lon <= -20:
        return True
    return False, highlandName

def highland2(lat, lon):
    # Highland 2
    highlandName = "Highland 2"
    if lat >= 20 and lat <= 40 and lon >= -40 and lon <= -20:
        return True
    return False, highlandName

def highland3(lat, lon):
    # Highland 3
    highlandName = "Highland 3"
    if lat >= -20 and lat <= 20 and lon >= -60 and lon <= -40:
        return True
    return False, highlandName

def in_highland(lat, lon):
    for highland in highland_regions:
        inHighland, highlandName = highland(lat, lon)
        if highland(lat, lon):
            return True, highlandName
    return False, None

# key values pairs with highland names as keys and highland functions as values
highland_regions = {
    "Highland 1": highland1,
    "Highland 2": highland2,
    "Highland 3": highland3
}
