def download_PCNet(net_id):
    import ndex2.client
    ndex = ndex2.client.Ndex2("http://public.ndexbio.org")
    summary = ndex.get_network_summary(net_id)
    version = summary['version']
    print(summary)
    cx_json = ndex.get_network_as_cx_stream(net_id).json()
    
    nodes = {}
    for aspect in cx_json:
        if not 'nodes' in aspect: continue
        for v in aspect['nodes']:
            nodes[v['@id']] = v['n']

    with open(f"data/PCNet_v{version}.tsv", "w") as oup:
        print("from\tto", file=oup)
        for aspect in cx_json:
            if not 'edges' in aspect: continue
            for e in aspect['edges']:
                print(nodes[e['s']] + "\t" + nodes[e['t']], file=oup)
    
    
if __name__ == "__main__":
    NET_ID = "4de852d9-9908-11e9-bcaf-0ac135e8bacf"
    download_PCNet(NET_ID)
