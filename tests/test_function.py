from txnova.function import _parse_hits


def test_parse_hits_keeps_confident_sdr() -> None:
    payload = {
        "results": [
            {
                "db": "afdb50",
                "alignments": [
                    [
                        {
                            "target": "AF-X-F1 Polyprenol dehydrogenase",
                            "prob": 1.0,
                            "eval": 1e-40,
                            "seqId": 50.6,
                            "taxName": "Mus musculus",
                        },
                        {
                            "target": "AF-Y-F1 junk",
                            "prob": 0.1,
                            "eval": 1.0,
                            "seqId": 10.0,
                        },
                    ]
                ],
            }
        ]
    }
    hits = _parse_hits(payload)
    assert len(hits) == 1
    assert hits[0]["description"] == "Polyprenol dehydrogenase"
    assert hits[0]["db"] == "afdb50"
