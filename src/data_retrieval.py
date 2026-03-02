"""
    methods for downloading files from S3 database

"""

import os
import boto3


def list_files(s3, bucket, folder): #, word_list):
    """
    lists files available which contain the words in word_list

    Args:
        s3 (boto3.client): client for s3
        bucket (str): bucket name
        folder (str): bucket folder
        word_list (list): list of relevant words (e.g. "OH2025_038")

    Returns:
        list: paths to relevant files
    """
    response = s3.list_objects_v2(Bucket=bucket, Prefix=folder)
    files = [content['Key'] for content in response['Contents']]
    # files = [file for file in files if any(word in file for word in word_list)]
    return files


def download_latest_from_s3(bucket, folder, target_path):
    """
    Downloads the most recently modified file in s3://bucket/folder/
    to the given local target_path.
    Uses AWS credentials from environment variables.
    """

    s3 = boto3.client(
        "s3",
        aws_access_key_id=os.getenv("AWS_ACCESS_KEY_ID"),
        aws_secret_access_key=os.getenv("AWS_SECRET_ACCESS_KEY"),
        aws_session_token=os.getenv("AWS_SESSION_TOKEN"),
    )

    # List all objects under the folder prefix
    response = s3.list_objects_v2(Bucket=bucket, Prefix=folder)

    if "Contents" not in response or len(response["Contents"]) == 0:
        raise FileNotFoundError(f"No files found in s3://{bucket}/{folder}")

    # Pick the newest file by LastModified timestamp
    latest_obj = max(response["Contents"], key=lambda x: x["LastModified"])
    key = latest_obj["Key"]

    # Download it
    filename = os.path.basename(key)   # the S3 filename

    os.makedirs(target_path, exist_ok=True)

    local_path = os.path.join(target_path, filename)

    s3.download_file(bucket, key, local_path)
    # s3.download_file(bucket, key, target_path)
    print('File downloaded')
    return key


def download_all_files_from_s3(bucket, folder, target_path): #, file_list):
    s3 = boto3.client(
        "s3",
        aws_access_key_id=os.getenv("AWS_ACCESS_KEY_ID"),
        aws_secret_access_key=os.getenv("AWS_SECRET_ACCESS_KEY"),
        aws_session_token=os.getenv("AWS_SESSION_TOKEN"),
    )
    available_files = list_files(s3, bucket, folder) #, file_list)
    i = 1
    for file in available_files:
        
        print(f"{i} of {len(available_files)}")
        i += 1
        download_file_from_s3(s3, bucket, file, target_path + file.split('/')[-1])
        print('Done!')


def download_specific_file_from_s3(bucket, folder, target_path, metadata_key='', metadata_value=''): 
    s3 = boto3.client(
        "s3",
        aws_access_key_id=os.getenv("AWS_ACCESS_KEY_ID"),
        aws_secret_access_key=os.getenv("AWS_SECRET_ACCESS_KEY"),
        aws_session_token=os.getenv("AWS_SESSION_TOKEN"),
    )
    available_files = list_files(s3, bucket, folder) #, file_list)
    i = 1
    if metadata_key and metadata_value:
        filtered_files = []
        print(f"Filtering files by metadata: {metadata_key} contains '{metadata_value}'")
        
        for file in available_files:
            try:
                # Get the object's metadata
                response = s3.head_object(Bucket=bucket, Key=file)
                metadata = response.get('Metadata', {})
                
                # Check if the metadata key exists and contains the value
                if metadata_key in metadata and metadata_value in metadata[metadata_key]:
                    filtered_files.append(file)
            except Exception as e:
                print(f"Error getting metadata for {file}: {e}")
        
        available_files = filtered_files
        print(f"Found {len(available_files)} files matching metadata filter")
        number_of_available_files = len(available_files)
        file_num = ''
        for file in available_files:
            print(f"{i} of {number_of_available_files}")
            if number_of_available_files>1:
                file_num = str(i)
            split_file = file.split('/')
            download_file_from_s3(s3, bucket, file, target_path + split_file[-1]+file_num)
            i += 1
        print('Done!')


def update_database(bucket):
    paths = [{"bucket_folder": "adat/parquet/base",
              "local_folder": "./data/adat/base/"},
              {"bucket_folder": "adat/parquet/anmlSMP",
              "local_folder": "./data/adat/anmlSMP/"},
              {"bucket_folder": "lab_poller/success",
              "local_folder": "./data/plates/"},]
    for path in paths:
        download_latest_from_s3(bucket, path['bucket_folder'], path['local_folder'])
    print("Done!")

def download_file_from_s3(s3, bucket, file, local_path):
    s3.download_file(bucket, file, local_path)


def download_plate(plate_name,
              bucket="com.oncohost-ops-prod-proteomic-db-49c25b67",
              destination='./data/plates',
              ):
    download_specific_file_from_s3(bucket,
                                   "lab_poller/success",
                                   destination,
                                   'plate_name',
                                   plate_name)
    
def download_adats(adat_stage="base", 
                   bucket="com.oncohost-ops-prod-proteomic-db-49c25b67",
                   destination='./data/plates',):
    baseline_path = "adat/parquet/"
    # stage_dict = {"base": "adat/parquet/",
    #               "anmlSMP": ""}
    
    download_latest_from_s3(bucket, baseline_path + adat_stage, destination)
    pass

if __name__ == '__main__':
    # update_database("com.oncohost-ops-prod-proteomic-db-49c25b67")
    # download_all_files_from_s3("com.oncohost-ops-prod-proteomic-db-49c25b67", "lab_poller/success", './data/plates/')
    # download_specific_file_from_s3("com.oncohost-ops-prod-proteomic-db-49c25b67", 
    #                                "lab_poller/success", 
    #                                './data/plates/', 
    #                                'plate_name',
    #                                'OH2025_047')


    plate_list = []
    # import os
    # print(os.listdir())
    plate_list += [i for i in range(49, 54, 1)]
    plate_list = ['OH2025_0' + str(pl) for pl in plate_list]
    plate_list = ['OH2025_055']
    for pl in plate_list:
        download_specific_file_from_s3("com.oncohost-ops-prod-proteomic-db-49c25b67", 
                                   "lab_poller/success", 
                                   './comparing izhar roee/data/', 
                                   'plate_name',
                                   pl)


def latest_file(folder_path):

    # Get all files in the folder
    files = os.listdir(folder_path)
    # Filter out directories, keep only files
    files = [folder_path + f for f in files if os.path.isfile(folder_path + f)]

    # Get the newest file by modification time
    newest_file = max(files, key=os.path.getmtime)

    # Get just the filename (without path)
    newest_filename = os.path.basename(newest_file)

    return newest_filename
    

    # download_latest_from_s3("com.oncohost-ops-prod-proteomic-db-49c25b67",
    #                          'adat/parquet/base', r'./data/adat/parquet/base/')
    # download_files_from_s3("com.oncohost-ops-prod-proteomic-db-49c25b67",
    #                         'adat/parquet/base',
    #                           './data/adat/parquet/base/mdacc_streck/',
    #                             ['OH2025_038', 'OH2025_039', 'OH2025_041', 'OH2025_043',
    #    'OH2025_047'])
