import pytest
from obspy.clients.fdsn import Client
from obspy.core.inventory.util import Latitude, Longitude
from pycheron.util.getLatLon import get_latlon


@pytest.fixture
def client_asset():
    return Client("IRIS")


def test_get_latlon_types(file_assets, client_asset):
    st = file_assets["7A_CABN_stream"]
    tr = st[0]
    lat, lon = get_latlon(client_asset, tr.stats.network, tr.stats.station)
    assert isinstance(lat, Latitude)
    assert isinstance(lon, Longitude)


def test_get_latlon_no_start_or_end(file_assets, client_asset):
    st = file_assets["7A_CABN_stream"]
    tr = st[0]
    lat, lon = get_latlon(client_asset, tr.stats.network, tr.stats.station)
    assert lat == 38.719898
    assert lon == -79.4412


def test_get_latlon_with_start_and_end(file_assets, client_asset):
    st = file_assets["7A_CABN_stream"]
    tr = st[0]
    lat, lon = get_latlon(
        client_asset,
        tr.stats.network,
        tr.stats.station,
        tr.stats.starttime.isoformat(),
        tr.stats.endtime.isoformat(),
    )
    assert lat == 38.719898
    assert lon == -79.4412


def test_get_latlon_failure(client_asset):
    lat, lon = get_latlon(
        client_asset,
        "BW",
        "BGLD",
    )
    assert lat is None
    assert lon is None
