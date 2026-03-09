from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

def protect_amine(input_smiles: str):
    """
    入力されたSMILESのアミノ基をアセチル基で保護し、
    結果のSMILESとハイライトされたSVG画像を返す。
    """
    # 1. SMILESから分子オブジェクトを作成（Boltz互換の立体情報等も保持）
    mol = Chem.MolFromSmiles(input_smiles)
    if mol is None:
        raise ValueError("無効なSMILES文字列です。")

    # 2. 反応ルールの定義 (SMARTS)
    # ターゲット: 1級または2級アミン [N;H2,H1:1]
    # 試薬: アセチルクロリド [CH3:2][C:3](=[O:4])[Cl]
    # 生成物: アミド結合 [N:1][C:3](=[O:4])[CH3:2]
    rxn_smarts = '[N;H2,H1:1].[CH3:2][C:3](=[O:4])[Cl]>>[N:1][C:3](=[O:4])[CH3:2]'
    rxn = AllChem.ReactionFromSmarts(rxn_smarts)
    
    # 試薬の準備
    reagent = Chem.MolFromSmiles('CC(=O)Cl')

    # 3. 反応の実行
    products = rxn.RunReactants((mol, reagent))
    
    if not products:
        return {"error": "この分子には反応するアミノ基がありません。"}

    # 最初の生成物を取得し、化学的に正しい状態にする（サニタイズ）
    product_mol = products[0][0]
    Chem.SanitizeMol(product_mol)
    
    # 立体化学などの情報をすべて保持したままSMILESを出力
    product_smiles = Chem.MolToSmiles(product_mol, isomericSmiles=True)

    # 4. 視覚化：新しく付加された保護基（アセチル基部分）を検索してハイライト
    # ※ハイライト用の部分構造SMARTS（N-C(=O)-C）
    highlight_pattern = Chem.MolFromSmarts('NC(=O)C')
    match_indices = product_mol.GetSubstructMatch(highlight_pattern)

    # SVG画像の生成 (Web表示用)
    drawer = rdMolDraw2D.MolDraw2DSVG(500, 300)
    drawer.DrawMolecule(product_mol, highlightAtoms=match_indices)
    drawer.FinishDrawing()
    svg_text = drawer.GetDrawingText()

    return {
        "original_smiles": input_smiles,
        "product_smiles": product_smiles,
        "svg_image": svg_text
    }
