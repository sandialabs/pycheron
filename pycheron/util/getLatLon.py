from obspy import UTCDateTime


def get_latlon(client, network, station, start=None, end=None):
    """
    Obtains latitude and longitude station info from an obspy client.

    :param client: Obspy Client object
    :type client: obspy.clients.fdsn.Client
    :param network: Network name from snclq id of obspy trace
    :type network: str
    :param station: Station name from snclq id of obspy trace
    :type station: str
    :param start: datetime string representation of starting time from collected obspy trace data
    :type start: str
    :param end: datetime string representation of ending time from collected obspy trace data
    :type end: str

    :return: latitude and longitude values from obspy station
    :rtype: (obspy.core.inventory.util.Latitude, obspy.core.inventory.util.Longitude)
    """
    try:
        if start is None or end is None:
            inv = client.get_stations(
                network=network, station=station, level="response"
            )
        else:
            starttime, endtime = UTCDateTime(start), UTCDateTime(end)
            inv = client.get_stations(
                network=network,
                station=station,
                starttime=starttime,
                endtime=endtime,
                level="response",
            )
    except Exception as e:
        print(f"Unable to generate lat or lon values: {e}")
        return None, None
    # lat/lon information is found in station. Can not use get_coordinates() because this requires "channel" info
    # which isn't always provided when setting a database metric
    lat = inv._networks[0]._stations[0]._latitude
    lon = inv._networks[0]._stations[0]._longitude
    return lat, lon
