import os
import requests
import json
import argparse
import pandas as pd
from typing import Tuple
import io
import numpy as np


class CVDSpider(object):
    headers = {
        'accept': '*/*',
        'accept-language': 'zh-CN,zh;q=0.9',
        'origin': 'https://cvd.hugeamp.org',
        'priority': 'u=1, i',
        'sec-ch-ua': '"Google Chrome";v="129", "Not=A?Brand";v="8", "Chromium";v="129"',
        'sec-ch-ua-mobile': '?0',
        'sec-ch-ua-platform': '"Windows"',
        'sec-fetch-dest': 'empty',
        'sec-fetch-mode': 'cors',
        'sec-fetch-site': 'same-site',
        'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/129.0.0.0 Safari/537.36',
        'Host': 'bioindex.hugeamp.org',
        'Connection': 'keep-alive'
    }
    Session = requests.Session()
    Session.headers = headers

    def __init__(self):
        self.rs_url = "https://bioindex.hugeamp.org/api/bio/varIdLookup/{}"
        self.var_url = "https://bioindex.hugeamp.org/api/bio/query/variant?q={}"
        self.token_url = "https://bioindex.hugeamp.org/api/bio/query/variant-dataset-associations?q={}"
        self.cont_url = "https://bioindex.hugeamp.org/api/bio/cont?token={}"
        self.data_set_url = "https://bioindex.hugeamp.org/api/portal/datasets?q=a2f"

        self.dataset = self.get_dataset()

    def get_dataset(self) -> pd.DataFrame:
        dataset_resp = self.Session.get(self.data_set_url)
        dataset_json = dataset_resp.json()
        dataset_data = dataset_json['data']

        # 使用 StringIO 包装 JSON 字符串
        data_str = json.dumps(dataset_data)
        dataset_df = pd.read_json(io.StringIO(data_str))

        return dataset_df

    def get_varid(self, rsid: str) -> Tuple[str, str]:
        varid_resp = self.Session.get(self.rs_url.format(rsid))
        varid_json = varid_resp.json()
        varid = varid_json['data']["varid"]
        varid_encode = varid.replace(":", "%3A")
        return varid, varid_encode

    def get_var_info(self, varid: str) -> str:
        var_resp = self.Session.get(self.var_url.format(varid))
        var_json = var_resp.json()
        var_nearest_gene: list = var_json['data'][0]["nearest"]
        return var_nearest_gene

    def get_token(self, varid: str) -> str | pd.DataFrame:
        token_resp = self.Session.get(self.token_url.format(varid))
        token_json = token_resp.json()
        token: str = token_json["continuation"]
        count = token_json["count"]
        if (token is None) or (count > 10):
            cont_data: list = token_json["data"]

            # 使用 StringIO 包装 JSON 字符串
            data_str = json.dumps(cont_data)
            data_df = pd.read_json(io.StringIO(data_str))

            return data_df
        else:
            return token

    def get_cont(self, token: str) -> pd.DataFrame:
        cont_resp = self.Session.get(self.cont_url.format(token))
        cont_json = cont_resp.json()
        cont_data: list = cont_json["data"]

        # 使用 StringIO 包装 JSON 字符串
        data_str = json.dumps(cont_data)
        data_df = pd.read_json(io.StringIO(data_str))

        return data_df

    def get_maxpValue(self, data_df: pd.DataFrame, phenotype: tuple) -> Tuple[float, str]:
        # print(phenotype)
        data_df = data_df[data_df['phenotype'].isin(phenotype)]
        data_df = data_df.sort_values(by='pValue', ascending=True)
        max_pv = data_df.iloc[0]['pValue']
        dataset = data_df.iloc[0]['dataset']
        return max_pv, dataset

    def get_PMID(self, dataset: str) -> str:
        pmid = self.dataset[self.dataset['name'] == dataset]['pmid'].values
        if all(pd.isnull(pmid)):
            pmid = dataset
        else:
            pmid = int(pmid[0])
        return str(pmid)


def main(rslist_file: str, phenotype: str, output_file: str, failed_output_file: str):
    CVDSpider_obj = CVDSpider()

    # 读取rslist文件中的SNP rsid
    with open(rslist_file, 'r') as f:
        rsids = [line.strip() for line in f.readlines()]

    # 检查输出文件是否存在，读取已处理的SNP
    if os.path.exists(output_file):
        processed_df = pd.read_csv(output_file)
        processed_rsids = set(processed_df['rsid'].tolist())
        print(f"发现已有输出文件，已处理的 rsid 数量为: {len(processed_rsids)}")
    else:
        # 如果输出文件不存在，创建一个空文件
        processed_rsids = set()
        with open(output_file, 'w') as f:
            # 写入CSV文件的表头
            f.write('rsid,varid,gene,pvalue,pmid\n')

    # 检查是否已有失败的SNP文件，读取已记录的失败SNP
    if os.path.exists(failed_output_file):
        failed_df = pd.read_csv(failed_output_file)
        failed_rsids = set(failed_df['rsid'].tolist())
        print(f"发现已有失败的 rsid 数量为: {len(failed_rsids)}")
    else:
        failed_rsids = set()
        with open(failed_output_file, 'w') as f:
            # 写入CSV文件的表头
            f.write('rsid\n')

    # 对每个rsid进行处理
    for i, rsid in enumerate(rsids, 1):
        if rsid in processed_rsids:
            print(f"跳过已处理的 rsid: {rsid}")
            continue
        if rsid in failed_rsids:
            print(f"跳过已记录为失败的 rsid: {rsid}")
            continue

        print(f"正在处理第 {i}/{len(rsids)} 个 rsid: {rsid}")
        try:
            # 获取 varid 和 varid_encode
            varid, varid_encode = CVDSpider_obj.get_varid(rsid)
            # 获取最近的基因信息
            gene = CVDSpider_obj.get_var_info(varid_encode)
            # 获取 token 或者直接的数据
            token = CVDSpider_obj.get_token(varid_encode)
            if isinstance(token, pd.DataFrame):
                data_df = token
            else:
                data_df = CVDSpider_obj.get_cont(token)
            # 获取最小 p-value 和数据集
            max_pv, dataset = CVDSpider_obj.get_maxpValue(data_df, phenotype)
            # 获取 PMID
            pmid = CVDSpider_obj.get_PMID(dataset)
            # 将结果追加写入CSV文件
            with open(output_file, 'a') as f:
                f.write(f"{rsid},{varid},{gene},{max_pv},{pmid}\n")

        except Exception as e:
            print(f"处理 rsid {rsid} 时出错: {e}")
            # 发生错误时，记录到失败文件
            with open(failed_output_file, 'a') as f:
                f.write(f"{rsid}\n")

    print(f"所有SNP处理完成，结果已保存到 {output_file}")
    print(f"处理失败的SNP已记录到 {failed_output_file}")


if __name__ == "__main__":
    
    def parse_tuple(value):
        return tuple(value.split(','))
    
    parser = argparse.ArgumentParser(description="CVD爬虫工具，获取SNP相关信息")
    parser.add_argument('--rslist', required=True,
                        help="输入包含rsid列表的文件路径，每行一个rsid")
    parser.add_argument('--pheno', type=parse_tuple, required=True, help="表型名称，使用逗号分隔的字符串")
    parser.add_argument('--failed', required=True, help="错误的rsid记录文件路径")
    parser.add_argument('--output', required=True, help="输出文件路径（CSV格式）")

    args = parser.parse_args()

    main(args.rslist, args.pheno, args.output, args.failed)
