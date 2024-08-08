Python codes for Geiger's method for earthquake hypocentre location. 
Loc_eq.py - main script
Iterate.py - modules for Loc_eq.py
Loc_eq_f.py - modules and functions fot Iterate.py
stacoord.lst - file containing list of stations and their geolocation
Bulletins are named in the format YYYY-MM-DD.HH-MM which defines preliminary hypocenter time of the earthquake;
Structure of the file is:
magnitude hypocenter_depth_in_km id
epicenter_latitude epicenter_longitude
station_code+phase_name phase_pick_in_seconds