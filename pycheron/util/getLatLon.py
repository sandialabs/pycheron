from obspy import UTCDateTime
import toml
import os
from pathlib import Path

__all__ = ["get_latlon"]


def get_latlon(
    client,
    network,
    station,
    start=None,
    end=None,
    manual=False,
    wfdb_conn=None,
):
    """
    Obtains latitude and longitude station info from an ObsPy client.

    :param client: Obspy Client object
    :type client: obspy.clients.fdsn.Client
    :param network: Network name from snclq id of ObsPy Trace
    :type network: str
    :param station: Station name from snclq id of ObsPy Trace
    :type station: str
    :param start: datetime string representation of starting time from collected ObsPy Trace data
    :type start: str
    :param end: datetime string representation of ending time from collected ObsPy Trace data
    :type end: str
    :param manual: If True, returns 0.0 for lat lon values.
                   Used for when data isn't available on IRIS so endpoints aren't repeatedly hit/failed.
    :type manual: bool
    :param wfdb_conn: database connection to query site table for lat/lon info
    :type wfdb_conn: WfdiscDB


    :return: latitude and longitude values from obspy station
    :rtype: (obspy.core.inventory.util.Latitude, obspy.core.inventory.util.Longitude)
    """

    # If database connection
    if wfdb_conn and not manual:
        # Try reading from the site table and getting the lat, lon
        try:
            query_result = wfdb_conn.read_site_table_location(wfdb_conn.table_names["site"], station, wfdb_conn.timebox)
            row = query_result.next()
            lat, lon = row["lat"], row["lon"]
        # If exception, set lat, lon to 0
        except Exception as e:
            print(f"Unable to extract table info. Assigning lat/lon as 0,0: {e}")
            lat, lon = 0.0, 0.0
        return lat, lon
    # If manual load of lat, lon information
    elif manual:
        # Get lat lon config file
        file_path = Path(os.path.realpath(__file__))
        conf_path = file_path.parent.parent.parent
        conf_file = Path(conf_path, "latlon_config.toml")
        conf = toml.load(conf_file)
        # Try to get lat, lon for specific net, sta
        try:
            lat = conf[network][station]["lat"]
            lon = conf[network][station]["lon"]
        # If excpetion, set lat, lon to 0
        except Exception as e:
            print(f"Unable to find config for network {network} and station {station}. Assigning lat/lon as 0,0: {e}")
            lat, lon = 0.0, 0.0
        return lat, lon
    # If not database, or manual, get the metadata from IRIS
    else:
        # Try to get response information is start/end is None -- probably don't want to allow this actually
        try:
            if start is None or end is None:
                inv = client.get_stations(network=network, station=station, level="response")
            # Get response information so we can extract lat, lon values
            else:
                starttime, endtime = UTCDateTime(start), UTCDateTime(end)
                inv = client.get_stations(
                    network=network,
                    station=station,
                    starttime=starttime,
                    endtime=endtime,
                    level="response",
                )
        # If exception note unable to get lat, lon values
        except Exception as e:
            print(f"Unable to generate lat or lon values: {e}")
            return None, None
        # lat/lon information is found in station. Can not use get_coordinates() because this requires "channel" info
        # which isn't always provided when setting a database metric
        # Get lat, lon values from inventory object
        lat = inv._networks[0]._stations[0]._latitude
        lon = inv._networks[0]._stations[0]._longitude
        return lat, lon
