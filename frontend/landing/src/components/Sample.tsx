import { FileText, ArrowRight, Database } from 'lucide-react';

export function Sample() {
  return (
    <div className="max-w-7xl mx-auto space-y-12 animate-in fade-in duration-700">
      {/* Header */}
      <div className="text-center space-y-4">
        <h1 className="text-5xl bg-gradient-to-r from-purple-400 to-purple-600 bg-clip-text text-transparent">
          Sample Inputs & Outputs
        </h1>
        <p className="text-xl text-gray-400">
          Example data formats for each pipeline
        </p>
      </div>

      {/* Generative Screening */}
      <div className="space-y-6">
        <div className="flex items-center gap-4">
          <div className="p-3 bg-gradient-to-br from-purple-500 to-pink-500 rounded-lg">
            <FileText className="w-6 h-6 text-white" />
          </div>
          <h2 className="text-3xl text-white">Generative Screening</h2>
        </div>
        
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          {/* Input Table */}
          <div className="bg-gradient-to-br from-purple-900/30 to-purple-950/20 border border-purple-500/30 rounded-xl overflow-hidden">
            <div className="bg-purple-500/20 px-4 py-3 border-b border-purple-500/30">
              <h3 className="text-white flex items-center gap-2">
                <Database className="w-4 h-4" />
                Input Columns
              </h3>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full">
                <thead className="bg-purple-900/30">
                  <tr>
                    <th className="px-4 py-3 text-left text-purple-300 text-sm">Column</th>
                    <th className="px-4 py-3 text-left text-purple-300 text-sm">Example</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-purple-500/10">
                  <tr className="hover:bg-purple-500/5">
                    <td className="px-4 py-3 text-gray-300">moleculechemblid</td>
                    <td className="px-4 py-3 text-gray-400">CHEMBL1</td>
                  </tr>
                  <tr className="hover:bg-purple-500/5">
                    <td className="px-4 py-3 text-gray-300">canonicalsmiles</td>
                    <td className="px-4 py-3 text-gray-400">CCCCC1CCCCC1O</td>
                  </tr>
                  <tr className="hover:bg-purple-500/5">
                    <td className="px-4 py-3 text-gray-300">MW</td>
                    <td className="px-4 py-3 text-gray-400">150.2</td>
                  </tr>
                  <tr className="hover:bg-purple-500/5">
                    <td className="px-4 py-3 text-gray-300">LogP</td>
                    <td className="px-4 py-3 text-gray-400">2.3</td>
                  </tr>
                  <tr className="hover:bg-purple-500/5">
                    <td className="px-4 py-3 text-gray-300">NumHDonors</td>
                    <td className="px-4 py-3 text-gray-400">1</td>
                  </tr>
                  <tr className="hover:bg-purple-500/5">
                    <td className="px-4 py-3 text-gray-300">NumHAcceptors</td>
                    <td className="px-4 py-3 text-gray-400">1</td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>

          {/* Output Table */}
          <div className="bg-gradient-to-br from-green-900/30 to-green-950/20 border border-green-500/30 rounded-xl overflow-hidden">
            <div className="bg-green-500/20 px-4 py-3 border-b border-green-500/30">
              <h3 className="text-white flex items-center gap-2">
                <ArrowRight className="w-4 h-4" />
                Output Columns
              </h3>
            </div>
            <div className="overflow-x-auto">
              <table className="w-full">
                <thead className="bg-green-900/30">
                  <tr>
                    <th className="px-4 py-3 text-left text-green-300 text-sm">Column</th>
                    <th className="px-4 py-3 text-left text-green-300 text-sm">Example</th>
                  </tr>
                </thead>
                <tbody className="divide-y divide-green-500/10">
                  <tr className="hover:bg-green-500/5">
                    <td className="px-4 py-3 text-gray-300">genmolid</td>
                    <td className="px-4 py-3 text-gray-400">GEN0001</td>
                  </tr>
                  <tr className="hover:bg-green-500/5">
                    <td className="px-4 py-3 text-gray-300">canonicalsmiles</td>
                    <td className="px-4 py-3 text-gray-400">CCNCCCC1CCCCC1O</td>
                  </tr>
                  <tr className="hover:bg-green-500/5">
                    <td className="px-4 py-3 text-gray-300">MW</td>
                    <td className="px-4 py-3 text-gray-400">179.3</td>
                  </tr>
                  <tr className="hover:bg-green-500/5">
                    <td className="px-4 py-3 text-gray-300">LogP</td>
                    <td className="px-4 py-3 text-gray-400">2.5</td>
                  </tr>
                  <tr className="hover:bg-green-500/5">
                    <td className="px-4 py-3 text-gray-300">pIC50</td>
                    <td className="px-4 py-3 text-gray-400">7.80</td>
                  </tr>
                  <tr className="hover:bg-green-500/5">
                    <td className="px-4 py-3 text-gray-300">activityclass</td>
                    <td className="px-4 py-3 text-gray-400">active</td>
                  </tr>
                </tbody>
              </table>
            </div>
          </div>
        </div>
      </div>

      {/* Virtual Screening */}
      <div className="space-y-6">
        <div className="flex items-center gap-4">
          <div className="p-3 bg-gradient-to-br from-purple-500 to-blue-500 rounded-lg">
            <FileText className="w-6 h-6 text-white" />
          </div>
          <h2 className="text-3xl text-white">Virtual Screening</h2>
        </div>
        
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="bg-gradient-to-br from-purple-900/30 to-purple-950/20 border border-purple-500/30 rounded-xl p-6">
            <h3 className="text-white mb-4 flex items-center gap-2">
              <Database className="w-4 h-4" />
              Input Example
            </h3>
            <div className="space-y-2 text-gray-300">
              <p><span className="text-purple-300">moleculechemblid:</span> 1</p>
              <p><span className="text-purple-300">canonicalsmiles:</span> CCOCOC1CCCCC1</p>
              <p><span className="text-purple-300">MW:</span> 164.2</p>
              <p><span className="text-purple-300">LogP:</span> 2.8</p>
            </div>
          </div>

          <div className="bg-gradient-to-br from-green-900/30 to-green-950/20 border border-green-500/30 rounded-xl p-6">
            <h3 className="text-white mb-4 flex items-center gap-2">
              <ArrowRight className="w-4 h-4" />
              Output Example
            </h3>
            <div className="space-y-2 text-gray-300">
              <p><span className="text-green-300">moleculechemblid:</span> 5</p>
              <p><span className="text-green-300">canonicalsmiles:</span> CCCCC1CCCCC1O</p>
              <p><span className="text-green-300">pIC50:</span> 8.30</p>
              <p><span className="text-green-300">activityclass:</span> active</p>
            </div>
          </div>
        </div>
      </div>

      {/* Predict Properties */}
      <div className="space-y-6">
        <div className="flex items-center gap-4">
          <div className="p-3 bg-gradient-to-br from-purple-500 to-indigo-500 rounded-lg">
            <FileText className="w-6 h-6 text-white" />
          </div>
          <h2 className="text-3xl text-white">Predict Properties</h2>
        </div>
        
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="bg-gradient-to-br from-purple-900/30 to-purple-950/20 border border-purple-500/30 rounded-xl overflow-hidden">
            <div className="bg-purple-500/20 px-4 py-3 border-b border-purple-500/30">
              <h3 className="text-white flex items-center gap-2">
                <Database className="w-4 h-4" />
                Input Columns
              </h3>
            </div>
            <div className="p-4 space-y-2 text-gray-300">
              <p><span className="text-purple-300">moleculechemblid:</span> idea1</p>
              <p><span className="text-purple-300">canonicalsmiles:</span> CCO</p>
              <p><span className="text-purple-300">MW:</span> 46.07</p>
              <p><span className="text-purple-300">LogP:</span> -0.3</p>
              <p><span className="text-purple-300">NumHDonors:</span> 1</p>
              <p><span className="text-purple-300">NumHAcceptors:</span> 1</p>
            </div>
          </div>

          <div className="bg-gradient-to-br from-green-900/30 to-green-950/20 border border-green-500/30 rounded-xl overflow-hidden">
            <div className="bg-green-500/20 px-4 py-3 border-b border-green-500/30">
              <h3 className="text-white flex items-center gap-2">
                <ArrowRight className="w-4 h-4" />
                Output Columns
              </h3>
            </div>
            <div className="p-4 space-y-2 text-gray-300">
              <p><span className="text-green-300">moleculechemblid:</span> idea1</p>
              <p><span className="text-green-300">canonicalsmiles:</span> CCO</p>
              <p><span className="text-green-300">MW:</span> 46.07</p>
              <p><span className="text-green-300">LogP:</span> -0.3</p>
              <p><span className="text-green-300">NumHDonors:</span> 1</p>
              <p><span className="text-green-300">NumHAcceptors:</span> 1</p>
              <p className="pt-2 border-t border-green-500/20">
                <span className="text-green-400">pIC50:</span> 6.20
              </p>
              <p><span className="text-green-400">activityclass:</span> active</p>
              <p><span className="text-green-400">solubility:</span> moderate</p>
            </div>
          </div>
        </div>
      </div>

      {/* Prioritize Candidates */}
      <div className="space-y-6">
        <div className="flex items-center gap-4">
          <div className="p-3 bg-gradient-to-br from-purple-500 to-violet-500 rounded-lg">
            <FileText className="w-6 h-6 text-white" />
          </div>
          <h2 className="text-3xl text-white">Prioritize Candidates</h2>
        </div>
        
        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <div className="bg-gradient-to-br from-purple-900/30 to-purple-950/20 border border-purple-500/30 rounded-xl overflow-hidden">
            <div className="bg-purple-500/20 px-4 py-3 border-b border-purple-500/30">
              <h3 className="text-white flex items-center gap-2">
                <Database className="w-4 h-4" />
                Input Columns
              </h3>
            </div>
            <div className="p-4 space-y-2 text-gray-300">
              <p><span className="text-purple-300">moleculechemblid:</span> idea1</p>
              <p><span className="text-purple-300">canonicalsmiles:</span> CCO</p>
              <p><span className="text-purple-300">pIC50:</span> 6.20</p>
              <p><span className="text-purple-300">activityclass:</span> active</p>
              <p><span className="text-purple-300">priorityrank:</span> 3</p>
            </div>
          </div>

          <div className="bg-gradient-to-br from-green-900/30 to-green-950/20 border border-green-500/30 rounded-xl overflow-hidden">
            <div className="bg-green-500/20 px-4 py-3 border-b border-green-500/30">
              <h3 className="text-white flex items-center gap-2">
                <ArrowRight className="w-4 h-4" />
                Output Columns (Sorted by Priority)
              </h3>
            </div>
            <div className="p-4 space-y-2 text-gray-300">
              <p><span className="text-green-300">moleculechemblid:</span> idea5</p>
              <p><span className="text-green-300">canonicalsmiles:</span> CCOCCO</p>
              <p><span className="text-green-300">pIC50:</span> 8.11</p>
              <p><span className="text-green-300">activityclass:</span> active</p>
              <p className="flex items-center gap-2">
                <span className="text-green-400">priorityrank:</span> 
                <span className="px-2 py-1 bg-green-500/30 rounded text-green-300 text-sm">1</span>
              </p>
            </div>
          </div>
        </div>
      </div>

      {/* Info Box */}
      <div className="bg-gradient-to-r from-purple-900/40 to-purple-800/40 border border-purple-500/40 rounded-2xl p-8">
        <h3 className="text-2xl text-white mb-4">Data Format Guidelines</h3>
        <ul className="space-y-2 text-gray-300">
          <li className="flex items-start gap-2">
            <span className="text-purple-400 mt-1">•</span>
            <span>Upload data in CSV format with appropriate column headers</span>
          </li>
          <li className="flex items-start gap-2">
            <span className="text-purple-400 mt-1">•</span>
            <span>SMILES notation should be in canonical format</span>
          </li>
          <li className="flex items-start gap-2">
            <span className="text-purple-400 mt-1">•</span>
            <span>Output files will include all input columns plus AI-generated predictions</span>
          </li>
          <li className="flex items-start gap-2">
            <span className="text-purple-400 mt-1">•</span>
            <span>Each pipeline processes data differently - refer to specific examples above</span>
          </li>
        </ul>
      </div>
    </div>
  );
}
