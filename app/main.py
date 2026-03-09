from fastapi import FastAPI, Request, Form
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from app.chem_engine import protect_amine

app = FastAPI()
templates = Jinja2Templates(directory="templates")

@app.get("/", response_class=HTMLResponse)
async def read_root(request: Request):
    # 初期画面の表示
    return templates.TemplateResponse("index.html", {"request": request})

@app.post("/predict", response_class=HTMLResponse)
async def predict_reaction(request: Request, smiles: str = Form(...)):
    try:
        # chem_engineを呼び出して反応を予測
        result = protect_amine(smiles)
        
        if "error" in result:
            return templates.TemplateResponse("index.html", {"request": request, "error": result["error"], "smiles": smiles})
            
        return templates.TemplateResponse("index.html", {
            "request": request,
            "original_smiles": result["original_smiles"],
            "product_smiles": result["product_smiles"],
            "svg_image": result["svg_image"]
        })
    except Exception as e:
        return templates.TemplateResponse("index.html", {"request": request, "error": str(e), "smiles": smiles})
