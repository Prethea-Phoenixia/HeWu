def _uc_mi2kft(dist):
    """
    Converts distance in miles to distance in kilofeet
    dist: distance in miles
    """
    conversion_const = 5.28  # kilofeet per mile
    return dist * conversion_const


def _uc_km2kft(dist):
    """
    Converts distance in kilometers to distance in kilofeet
    dist: distance in kilometers
    """
    conversion_const = 3.28084  # kilofeet per kilometer
    return dist * conversion_const


def _uc_kft2km(dist):
    """
    Converts distance in kilofeet to distance in kilometers
    dist: distance in kilofeet
    """
    conversion_const = 1 / 3.28084  # kilometer per kilofeet
    return dist * conversion_const


def _uc_kft2mi(dist):
    """
    Converts distance in kilofeet to distance in kilometers
    dist: distance in kilofeet
    """
    conversion_const = 1 / 5.28  # miles per kilofeet
    return dist * conversion_const


def _uc_mi2ft(dist):
    """
    Converts distance in miles to distance in feet
    dist: distance in miles
    """
    conversion_const = 1 / 5280  # miles per feet
    return dist * conversion_const


def _uc_ft2mi(dist):
    """
    Converts distance in feet to distance in miles
    dist: distance in feet
    """
    conversion_const = 5280  # feet per mile
    return dist * conversion_const


def _uc_km2ft(dist):
    """
    Converts distance in kilometers to distance in feet
    dist: distance in kilometers
    """
    conversion_const = 3280.84  # feet per kilometer
    return dist * conversion_const


def _uc_ft2km(dist):
    """
    Converts distance in feet to distance in kilometers
    dist: distance in feet
    """
    conversion_const = 1 / 3280.84  # kilometers per foot
    return dist * conversion_const


def _uc_m2kft(dist):
    """
    Converts distance in meter to kilofeet
    dist: distance in meter
    """
    conversion_const = 3.2808e-3  # kilo feet per meter
    return dist * conversion_const


def _uc_psi2pa(pres):
    """
    Converts presure from psi to pa
    pres: pressure in psi
    """
    conversion_const = 6894.76  # pa per psi
    return pres * conversion_const


def _uc_pa2psi(pres):
    """
    Converts presure from pa to psi
    pres: pressure in pa
    """
    conversion_const = 6894.76  # pa per psi
    return pres * conversion_const


def _uc_m2ft(dist):
    """
    Converts distance in feet to distance in meter
    dist: distance in feet
    """
    conversion_const = 3.28084  # feet per meter
    return dist * conversion_const


def _uc_ft2m(dist):
    """
    Converts distance in feet to distance in meter
    dist: distance in feet
    """
    conversion_const = 1 / 3.28084  # meter per foot
    return dist * conversion_const
