import json, os, logging

def download_PCNet(config):
    import ndex2.client
    ndex = ndex2.client.Ndex2(config['url'])
    summary = ndex.get_network_summary(config['id'])
    version = summary['version']
    cx_json = ndex.get_network_as_cx_stream(config['id']).json()
    
    nodes = {}
    for aspect in cx_json:
        if not 'nodes' in aspect: continue
        for v in aspect['nodes']:
            nodes[v['@id']] = v['n']

    cur_dir = os.path.dirname(__file__)
    DOWNLOAD_PATH = f"{cur_dir}/data/PCNet_v{version}.tsv"
    with open(DOWNLOAD_PATH, "w") as oup:
        print("from\tto", file=oup)
        for aspect in cx_json:
            if not 'edges' in aspect: continue
            for e in aspect['edges']:
                print(nodes[e['s']] + "\t" + nodes[e['t']], file=oup)

def download_DisGeNET(config):
    from google_drive_downloader import GoogleDriveDownloader as gdd

    cur_dir = os.path.dirname(__file__)

    for f in config:
        dest_path = f"{cur_dir}/data/{f['filename']}"
        gdd.download_file_from_google_drive(file_id=f['id'],
                                        dest_path=dest_path,
                                        unzip=False)


def main():
    FORMAT = '%(asctime)-15s %(name)s %(levelname)s %(message)s'    
    logging.basicConfig(level=logging.INFO, format=FORMAT)    
    
    CONFIG_PATH =  os.path.join(os.path.dirname(__file__), 
                                "data/download_config.json")
    with open(CONFIG_PATH, "r") as json_file:
        configure = json.load(json_file)
    
    cur_dir = os.path.dirname(__file__)
    db_path = os.path.join(cur_dir, "data")
    db_path = os.path.abspath(db_path)
    if not os.path.exists(db_path):
        logging.info(f"Creating {db_path} to store outputs.")
        os.makedirs(db_path)    
    else:
        logging.info(f"Skipping creating a directory to store files.")
    
    logging.info("downloading PCNet...")
    download_PCNet(configure["PCNet"])
    
    logging.info("downloading DisGeNET...")    
    download_DisGeNET(configure["DisGeNET"])
    
    logging.info(f"Files have been downloaded to {db_path}.")