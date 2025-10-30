import { FileText, ArrowRight, Database, Home, Brain, Download, Lightbulb, Filter, Target, ListOrdered } from 'lucide-react';
import { Link } from 'react-router-dom';

export function SamplePage() {
  return (
    <div className="min-h-screen bg-[#0a0a0f] text-white">
      {/* Gradient Backgrounds */}
      <div className="fixed inset-0 overflow-hidden pointer-events-none">
        <div className="absolute top-0 right-0 w-[600px] h-[600px] bg-indigo-600/20 rounded-full blur-[120px]" />
        <div className="absolute bottom-0 left-0 w-[500px] h-[500px] bg-purple-600/15 rounded-full blur-[100px]" />
      </div>

      <div className="relative z-10">
        {/* Header */}
        <header className="px-8 py-6 flex items-center justify-between backdrop-blur-sm bg-white/5 border-b border-white/10 sticky top-0 z-50">
          <Link to="/" className="flex items-center gap-3 hover:opacity-80 transition-opacity group">
            <div className="p-2 bg-indigo-500/20 rounded-lg group-hover:bg-indigo-500/30 transition-colors">
              <Home className="w-5 h-5 text-indigo-400" />
            </div>
            <span className="text-white/80">Back to Home</span>
          </Link>
          
          <div className="flex items-center gap-3">
            <div className="p-2 bg-gradient-to-br from-indigo-500 to-indigo-700 rounded-lg">
              <Brain className="w-6 h-6 text-white" />
            </div>
            <div className="flex flex-col">
              <span className="text-xl bg-gradient-to-r from-indigo-400 to-indigo-600 bg-clip-text text-transparent">
                AI-DRIVEN DDD
              </span>
            </div>
          </div>
          
          <div className="flex items-center gap-4">
            <Link to="/overview" className="px-6 py-2 rounded-full bg-white/5 text-white/70 hover:bg-white/10 transition-all border border-white/10">
              Overview
            </Link>
            <Link to="/uses" className="px-6 py-2 rounded-full bg-white/5 text-white/70 hover:bg-white/10 transition-all border border-white/10">
              Uses
            </Link>
          </div>
        </header>

        {/* Content */}
        <main className="container mx-auto px-8 py-16">
          <div className="max-w-7xl mx-auto space-y-16">
            {/* Header */}
            <div className="text-center space-y-6">
              <div className="inline-flex items-center gap-2 px-4 py-2 rounded-full bg-indigo-500/10 border border-indigo-500/20 mb-4">
                <FileText className="w-4 h-4 text-indigo-400" />
                <span className="text-sm text-indigo-300">Sample Data</span>
              </div>
              
              <h1 className="text-6xl leading-tight">
                <span className="block text-white mb-2">Sample</span>
                <span className="bg-gradient-to-r from-indigo-400 via-purple-500 to-purple-600 bg-clip-text text-transparent">
                  Inputs & Outputs
                </span>
              </h1>
              
              <p className="text-xl text-gray-400 max-w-3xl mx-auto">
                Example data formats for each pipeline with real-world examples
              </p>
            </div>

            {/* Generative Screening */}
            <div className="space-y-8">
              <div className="flex items-center gap-4">
                <div className="p-4 bg-gradient-to-br from-purple-500 to-pink-500 rounded-2xl shadow-lg">
                  <Lightbulb className="w-8 h-8 text-white" />
                </div>
                <div>
                  <h2 className="text-4xl text-white">Generative Screening</h2>
                  <p className="text-gray-400">AI-generated molecules from seed compounds</p>
                </div>
              </div>
              
              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                {/* Input Table */}
                <div className="bg-gradient-to-br from-purple-950/40 to-purple-900/20 border border-purple-500/30 rounded-2xl overflow-hidden hover:shadow-xl hover:shadow-purple-500/20 transition-all">
                  <div className="bg-purple-500/20 px-6 py-4 border-b border-purple-500/30 flex items-center justify-between">
                    <h3 className="text-white flex items-center gap-2 text-lg">
                      <Database className="w-5 h-5" />
                      Input Columns
                    </h3>
                    <span className="px-3 py-1 bg-purple-400/20 rounded-full text-sm text-purple-300">CSV</span>
                  </div>
                  <div className="overflow-x-auto">
                    <table className="w-full">
                      <thead className="bg-purple-900/40">
                        <tr>
                          <th className="px-6 py-4 text-left text-purple-300">Column</th>
                          <th className="px-6 py-4 text-left text-purple-300">Example</th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-purple-500/10">
                        <tr className="hover:bg-purple-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">moleculechemblid</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">CHEMBL1</td>
                        </tr>
                        <tr className="hover:bg-purple-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">canonicalsmiles</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">CCCCC1CCCCC1O</td>
                        </tr>
                        <tr className="hover:bg-purple-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">MW</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">150.2</td>
                        </tr>
                        <tr className="hover:bg-purple-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">LogP</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">2.3</td>
                        </tr>
                        <tr className="hover:bg-purple-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">NumHDonors</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">1</td>
                        </tr>
                        <tr className="hover:bg-purple-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">NumHAcceptors</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">1</td>
                        </tr>
                      </tbody>
                    </table>
                  </div>
                </div>

                {/* Output Table */}
                <div className="bg-gradient-to-br from-green-950/40 to-green-900/20 border border-green-500/30 rounded-2xl overflow-hidden hover:shadow-xl hover:shadow-green-500/20 transition-all">
                  <div className="bg-green-500/20 px-6 py-4 border-b border-green-500/30 flex items-center justify-between">
                    <h3 className="text-white flex items-center gap-2 text-lg">
                      <ArrowRight className="w-5 h-5" />
                      Output Columns
                    </h3>
                    <span className="px-3 py-1 bg-green-400/20 rounded-full text-sm text-green-300">Generated</span>
                  </div>
                  <div className="overflow-x-auto">
                    <table className="w-full">
                      <thead className="bg-green-900/40">
                        <tr>
                          <th className="px-6 py-4 text-left text-green-300">Column</th>
                          <th className="px-6 py-4 text-left text-green-300">Example</th>
                        </tr>
                      </thead>
                      <tbody className="divide-y divide-green-500/10">
                        <tr className="hover:bg-green-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">genmolid</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">GEN0001</td>
                        </tr>
                        <tr className="hover:bg-green-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">canonicalsmiles</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">CCNCCCC1CCCCC1O</td>
                        </tr>
                        <tr className="hover:bg-green-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">MW</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">179.3</td>
                        </tr>
                        <tr className="hover:bg-green-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">LogP</td>
                          <td className="px-6 py-4 text-gray-400 font-mono text-sm">2.5</td>
                        </tr>
                        <tr className="hover:bg-green-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">pIC50</td>
                          <td className="px-6 py-4 text-green-400 font-mono text-sm">7.80</td>
                        </tr>
                        <tr className="hover:bg-green-500/5 transition-colors">
                          <td className="px-6 py-4 text-gray-300">activityclass</td>
                          <td className="px-6 py-4">
                            <span className="px-3 py-1 bg-green-500/20 rounded-full text-sm text-green-300">active</span>
                          </td>
                        </tr>
                      </tbody>
                    </table>
                  </div>
                </div>
              </div>
            </div>

            {/* Virtual Screening */}
            <div className="space-y-8">
              <div className="flex items-center gap-4">
                <div className="p-4 bg-gradient-to-br from-blue-500 to-cyan-500 rounded-2xl shadow-lg">
                  <Filter className="w-8 h-8 text-white" />
                </div>
                <div>
                  <h2 className="text-4xl text-white">Virtual Screening</h2>
                  <p className="text-gray-400">Filtered and ranked candidates</p>
                </div>
              </div>
              
              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                <div className="bg-gradient-to-br from-blue-950/40 to-blue-900/20 border border-blue-500/30 rounded-2xl p-8 hover:shadow-xl hover:shadow-blue-500/20 transition-all">
                  <h3 className="text-white mb-6 flex items-center gap-2 text-lg">
                    <Database className="w-5 h-5" />
                    Input Example
                  </h3>
                  <div className="space-y-3 text-gray-300">
                    <div className="flex justify-between py-2 border-b border-blue-500/10">
                      <span className="text-blue-300">moleculechemblid:</span>
                      <span className="font-mono text-sm">1</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-blue-500/10">
                      <span className="text-blue-300">canonicalsmiles:</span>
                      <span className="font-mono text-sm">CCOCOC1CCCCC1</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-blue-500/10">
                      <span className="text-blue-300">MW:</span>
                      <span className="font-mono text-sm">164.2</span>
                    </div>
                    <div className="flex justify-between py-2">
                      <span className="text-blue-300">LogP:</span>
                      <span className="font-mono text-sm">2.8</span>
                    </div>
                  </div>
                </div>

                <div className="bg-gradient-to-br from-green-950/40 to-green-900/20 border border-green-500/30 rounded-2xl p-8 hover:shadow-xl hover:shadow-green-500/20 transition-all">
                  <h3 className="text-white mb-6 flex items-center gap-2 text-lg">
                    <ArrowRight className="w-5 h-5" />
                    Output Example
                  </h3>
                  <div className="space-y-3 text-gray-300">
                    <div className="flex justify-between py-2 border-b border-green-500/10">
                      <span className="text-green-300">moleculechemblid:</span>
                      <span className="font-mono text-sm">5</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-green-500/10">
                      <span className="text-green-300">canonicalsmiles:</span>
                      <span className="font-mono text-sm">CCCCC1CCCCC1O</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-green-500/10">
                      <span className="text-green-300">pIC50:</span>
                      <span className="font-mono text-sm text-green-400">8.30</span>
                    </div>
                    <div className="flex justify-between py-2">
                      <span className="text-green-300">activityclass:</span>
                      <span className="px-3 py-1 bg-green-500/20 rounded-full text-sm text-green-300">active</span>
                    </div>
                  </div>
                </div>
              </div>
            </div>

            {/* Predict Properties */}
            <div className="space-y-8">
              <div className="flex items-center gap-4">
                <div className="p-4 bg-gradient-to-br from-indigo-500 to-purple-500 rounded-2xl shadow-lg">
                  <Target className="w-8 h-8 text-white" />
                </div>
                <div>
                  <h2 className="text-4xl text-white">Predict Properties</h2>
                  <p className="text-gray-400">ML-powered property predictions</p>
                </div>
              </div>
              
              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                <div className="bg-gradient-to-br from-indigo-950/40 to-indigo-900/20 border border-indigo-500/30 rounded-2xl overflow-hidden hover:shadow-xl hover:shadow-indigo-500/20 transition-all">
                  <div className="bg-indigo-500/20 px-6 py-4 border-b border-indigo-500/30">
                    <h3 className="text-white flex items-center gap-2 text-lg">
                      <Database className="w-5 h-5" />
                      Input Columns
                    </h3>
                  </div>
                  <div className="p-6 space-y-3 text-gray-300">
                    <div className="flex justify-between py-2 border-b border-indigo-500/10">
                      <span className="text-indigo-300">moleculechemblid:</span>
                      <span className="font-mono text-sm">idea1</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-indigo-500/10">
                      <span className="text-indigo-300">canonicalsmiles:</span>
                      <span className="font-mono text-sm">CCO</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-indigo-500/10">
                      <span className="text-indigo-300">MW:</span>
                      <span className="font-mono text-sm">46.07</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-indigo-500/10">
                      <span className="text-indigo-300">LogP:</span>
                      <span className="font-mono text-sm">-0.3</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-indigo-500/10">
                      <span className="text-indigo-300">NumHDonors:</span>
                      <span className="font-mono text-sm">1</span>
                    </div>
                    <div className="flex justify-between py-2">
                      <span className="text-indigo-300">NumHAcceptors:</span>
                      <span className="font-mono text-sm">1</span>
                    </div>
                  </div>
                </div>

                <div className="bg-gradient-to-br from-green-950/40 to-green-900/20 border border-green-500/30 rounded-2xl overflow-hidden hover:shadow-xl hover:shadow-green-500/20 transition-all">
                  <div className="bg-green-500/20 px-6 py-4 border-b border-green-500/30">
                    <h3 className="text-white flex items-center gap-2 text-lg">
                      <ArrowRight className="w-5 h-5" />
                      Output Columns
                    </h3>
                  </div>
                  <div className="p-6 space-y-3 text-gray-300">
                    <div className="flex justify-between py-2 border-b border-green-500/10 opacity-50">
                      <span className="text-green-300">moleculechemblid:</span>
                      <span className="font-mono text-sm">idea1</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-green-500/10 opacity-50">
                      <span className="text-green-300">canonicalsmiles:</span>
                      <span className="font-mono text-sm">CCO</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-green-500/10 opacity-50">
                      <span className="text-green-300">MW:</span>
                      <span className="font-mono text-sm">46.07</span>
                    </div>
                    <div className="pt-3 border-t-2 border-green-500/30">
                      <p className="text-xs text-green-400 mb-3">AI-Predicted Properties:</p>
                      <div className="flex justify-between py-2 border-b border-green-500/10">
                        <span className="text-green-400">pIC50:</span>
                        <span className="font-mono text-sm text-green-300">6.20</span>
                      </div>
                      <div className="flex justify-between py-2 border-b border-green-500/10">
                        <span className="text-green-400">activityclass:</span>
                        <span className="px-2 py-1 bg-green-500/20 rounded text-sm text-green-300">active</span>
                      </div>
                      <div className="flex justify-between py-2">
                        <span className="text-green-400">solubility:</span>
                        <span className="font-mono text-sm text-green-300">moderate</span>
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            </div>

            {/* Prioritize Candidates */}
            <div className="space-y-8">
              <div className="flex items-center gap-4">
                <div className="p-4 bg-gradient-to-br from-violet-500 to-fuchsia-500 rounded-2xl shadow-lg">
                  <ListOrdered className="w-8 h-8 text-white" />
                </div>
                <div>
                  <h2 className="text-4xl text-white">Prioritize Candidates</h2>
                  <p className="text-gray-400">Ranked by predicted performance</p>
                </div>
              </div>
              
              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
                <div className="bg-gradient-to-br from-violet-950/40 to-violet-900/20 border border-violet-500/30 rounded-2xl overflow-hidden hover:shadow-xl hover:shadow-violet-500/20 transition-all">
                  <div className="bg-violet-500/20 px-6 py-4 border-b border-violet-500/30">
                    <h3 className="text-white flex items-center gap-2 text-lg">
                      <Database className="w-5 h-5" />
                      Input Columns
                    </h3>
                  </div>
                  <div className="p-6 space-y-3 text-gray-300">
                    <div className="flex justify-between py-2 border-b border-violet-500/10">
                      <span className="text-violet-300">moleculechemblid:</span>
                      <span className="font-mono text-sm">idea1</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-violet-500/10">
                      <span className="text-violet-300">canonicalsmiles:</span>
                      <span className="font-mono text-sm">CCO</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-violet-500/10">
                      <span className="text-violet-300">pIC50:</span>
                      <span className="font-mono text-sm">6.20</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-violet-500/10">
                      <span className="text-violet-300">activityclass:</span>
                      <span className="px-2 py-1 bg-violet-500/20 rounded text-sm text-violet-300">active</span>
                    </div>
                    <div className="flex justify-between py-2">
                      <span className="text-violet-300">priorityrank:</span>
                      <span className="px-3 py-1 bg-violet-500/30 rounded text-violet-300">3</span>
                    </div>
                  </div>
                </div>

                <div className="bg-gradient-to-br from-green-950/40 to-green-900/20 border border-green-500/30 rounded-2xl overflow-hidden hover:shadow-xl hover:shadow-green-500/20 transition-all">
                  <div className="bg-green-500/20 px-6 py-4 border-b border-green-500/30 flex items-center justify-between">
                    <h3 className="text-white flex items-center gap-2 text-lg">
                      <ArrowRight className="w-5 h-5" />
                      Output (Sorted by Priority)
                    </h3>
                  </div>
                  <div className="p-6 space-y-3 text-gray-300">
                    <div className="flex justify-between py-2 border-b border-green-500/10">
                      <span className="text-green-300">moleculechemblid:</span>
                      <span className="font-mono text-sm">idea5</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-green-500/10">
                      <span className="text-green-300">canonicalsmiles:</span>
                      <span className="font-mono text-sm">CCOCCO</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-green-500/10">
                      <span className="text-green-300">pIC50:</span>
                      <span className="font-mono text-sm text-green-400">8.11</span>
                    </div>
                    <div className="flex justify-between py-2 border-b border-green-500/10">
                      <span className="text-green-300">activityclass:</span>
                      <span className="px-2 py-1 bg-green-500/20 rounded text-sm text-green-300">active</span>
                    </div>
                    <div className="flex justify-between py-2">
                      <span className="text-green-300">priorityrank:</span>
                      <div className="flex items-center gap-2">
                        <span className="px-3 py-1 bg-green-500/40 rounded text-green-200">1</span>
                        <span className="text-xs text-green-400">Top Priority</span>
                      </div>
                    </div>
                  </div>
                </div>
              </div>
            </div>

            {/* Guidelines */}
            <div className="bg-gradient-to-r from-purple-900/40 to-indigo-800/40 border border-purple-500/40 rounded-3xl p-10">
              <div className="flex items-center gap-3 mb-6">
                <Download className="w-6 h-6 text-purple-400" />
                <h3 className="text-3xl text-white">Data Format Guidelines</h3>
              </div>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                <div className="space-y-3">
                  <div className="flex items-start gap-3">
                    <div className="w-2 h-2 rounded-full bg-purple-400 mt-2" />
                    <span className="text-gray-300">Upload data in CSV format with appropriate column headers</span>
                  </div>
                  <div className="flex items-start gap-3">
                    <div className="w-2 h-2 rounded-full bg-purple-400 mt-2" />
                    <span className="text-gray-300">SMILES notation should be in canonical format</span>
                  </div>
                </div>
                <div className="space-y-3">
                  <div className="flex items-start gap-3">
                    <div className="w-2 h-2 rounded-full bg-purple-400 mt-2" />
                    <span className="text-gray-300">Output files include all input columns plus AI predictions</span>
                  </div>
                  <div className="flex items-start gap-3">
                    <div className="w-2 h-2 rounded-full bg-purple-400 mt-2" />
                    <span className="text-gray-300">Each pipeline processes data differently - refer to examples above</span>
                  </div>
                </div>
              </div>
            </div>
          </div>
        </main>
      </div>
    </div>
  );
}
